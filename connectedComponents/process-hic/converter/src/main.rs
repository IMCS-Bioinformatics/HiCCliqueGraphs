use clap::Parser;
use hic_lib::{read_data, write_structs_to_csv, ConverterError, UniversalInteraction};
use serde::Deserialize;
use std::{error::Error, fs};

#[derive(Parser, Debug)]
struct Args {
    /// Length of chromatin segments; used to calculate start and end of segments from midpoint
    #[arg(short, long, default_value_t = 5000)]
    segment_len: u32,

    /// P value of files to read
    #[arg(short, long)]
    file_pvalue: f32,

    /// P value by which to filter data; should be larger or equal to FILE_PVALUE
    #[arg(short, long)]
    pvalue: f32,

    #[arg(
        short,
        long,
        default_value = "..\\Source_data\\Normal_HiC(hg19,3DIV_legacy)\\"
    )]
    og_interaction_path: String,

    #[arg(
        short,
        long,
        default_value = "..\\Processed_data\\Normal_Hi-C\\Universal\\"
    )]
    universal_data_path: String,
}

// Interaction contains the main fields from original Hi-C interaction data.
#[derive(Debug, Deserialize)]
struct Interaction {
    chr: String,
    bin1: u32,
    bin2: u32,
    significance: f32,
}

// Returns interaction transformed into universal format
fn make_universal_record(rec: &Interaction, segm_len: u32) -> UniversalInteraction {
    let universal_record: UniversalInteraction = UniversalInteraction {
        chr: rec.chr.clone(),
        bin1_start: rec.bin1.checked_sub(segm_len / 2).unwrap_or(0),
        bin1_end: rec.bin1 + segm_len / 2,
        bin2_start: rec.bin2.checked_sub(segm_len / 2).unwrap_or(0),
        bin2_end: rec.bin2 + segm_len / 2,
        p_value: rec.significance,
    };
    universal_record
}

// Creates a sorted vector of universal format interaction data, filtering out non-significant interactions
fn transform_data_to_universal(
    data: &Vec<Interaction>,
    pvalue: f32,
    segm_len: u32,
) -> Vec<UniversalInteraction> {
    let mut universal_data = Vec::new();
    for record in data {
        if record.significance >= pvalue {
            universal_data.push(make_universal_record(&record, segm_len));
        }
    }
    universal_data.sort_by_key(|item| (item.chr.clone(), item.bin1_start, item.bin2_start));
    universal_data
}

fn create_output_filename(filename: &str) -> Result<String, Box<dyn Error>> {
    let directories = filename.split("\\").collect::<Vec<&str>>();
    let tissue_name: &str = directories
        .get(directories.len() - 2)
        .ok_or("could not get tissue name")?;
    let tissue_id: &str = directories
        .last()
        .ok_or("could not get last directory")?
        .split("_")
        .collect::<Vec<&str>>()
        .get(1)
        .ok_or("could not get tissue id")?;
    let output_filename = format!("{}_{}.csv", tissue_id, tissue_name);
    Ok(output_filename)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    // If PVALUE is smaller than FILE_PVALUE, the result will only include interactions with p-value >= FILE_PVALUE
    // but folder name will contain PVALUE, which is misleading
    if args.pvalue < args.file_pvalue {
        panic!("P value should not be smaller than file p value");
    }
    if !std::path::Path::new(&args.og_interaction_path).exists() {
        return Err(Box::new(ConverterError::FileError(format!(
            "input directory '{}' does not exist",
            &args.og_interaction_path
        ))));
    }
    let paths = fs::read_dir(&args.og_interaction_path)?;
    let output_directory = format!("{0}\\Pvalue{1}\\", args.universal_data_path, args.pvalue)
        .to_string()
        .to_owned();
    fs::create_dir_all(&output_directory)?;
    for path in paths {
        let foldername = path?.path();
        let file_paths = fs::read_dir(foldername)?;
        for file_entry in file_paths {
            let file_path = file_entry?.path();
            let mut og_data: Vec<Interaction> = Vec::new();
            let filename = file_path.to_str().ok_or(format!(
                "Problem processing file path; check if path {} is valid unicode",
                file_path.display()
            ))?;
            if !filename.ends_with(&format!("-{0}.zip", args.file_pvalue).to_string()) {
                continue;
            }
            read_data(&mut og_data, &filename, b'\t', None)?;
            let universal_data =
                transform_data_to_universal(&og_data, args.pvalue, args.segment_len);
            write_structs_to_csv(
                &universal_data,
                &output_directory,
                &create_output_filename(filename)?,
            )?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transform_data_to_universal() {
        let e7 = 10_u32.pow(7);
        let mut data = Vec::new();
        // Generate some data
        for i in 1..7 {
            for j in (1..7).rev() {
                let interaction: Interaction = Interaction {
                    chr: "chr".to_string() + &i.to_string(),
                    bin1: (10 + j - i) * e7,
                    bin2: (10 + j - i) * e7 + 20000,
                    significance: (9 - j) as f32,
                };
                data.push(interaction);
            }
        }
        // Add some data by hand
        data.push(Interaction {
            chr: "chrX".to_string(),
            bin1: 10 * e7,
            bin2: 13 * e7,
            significance: 7.0,
        });
        data.push(Interaction {
            chr: "chr3".to_string(),
            bin1: 10 * e7,
            bin2: 13 * e7,
            significance: 7.0,
        });
        data.push(Interaction {
            chr: "chr12".to_string(),
            bin1: 8 * e7,
            bin2: 9 * e7,
            significance: 5.0,
        });
        // Create expected result
        let mut expected_res = Vec::new();
        for i in 1..7 {
            // Only generate data with pvalue above threshold
            for j in 1..5 {
                let interaction: UniversalInteraction = UniversalInteraction {
                    chr: "chr".to_string() + &i.to_string(),
                    bin1_start: (10 + j - i) * e7 - 2500,
                    bin1_end: (10 + j - i) * e7 + 2500,
                    bin2_start: (10 + j - i) * e7 + 20000 - 2500,
                    bin2_end: (10 + j - i) * e7 + 20000 + 2500,
                    p_value: (9 - j) as f32,
                };
                expected_res.push(interaction);
                if i == 3 && j == 3 {
                    expected_res.push(UniversalInteraction {
                        chr: "chr3".to_string(),
                        bin1_start: 10 * e7 - 2500,
                        bin1_end: 10 * e7 + 2500,
                        bin2_start: 13 * e7 - 2500,
                        bin2_end: 13 * e7 + 2500,
                        p_value: 7.0,
                    });
                }
            }
            if i == 1 {
                expected_res.push(UniversalInteraction {
                    chr: "chr12".to_string(),
                    bin1_start: 8 * e7 - 2500,
                    bin1_end: 8 * e7 + 2500,
                    bin2_start: 9 * e7 - 2500,
                    bin2_end: 9 * e7 + 2500,
                    p_value: 5.0,
                });
            }
        }
        expected_res.push(UniversalInteraction {
            chr: "chrX".to_string(),
            bin1_start: 10 * e7 - 2500,
            bin1_end: 10 * e7 + 2500,
            bin2_start: 13 * e7 - 2500,
            bin2_end: 13 * e7 + 2500,
            p_value: 7.0,
        });

        let res = transform_data_to_universal(&data, 5.0, 5000);
        // Results should be the same and in the same order
        // assert_eq!(res, expected_res);
        for (idx, record) in res.iter().enumerate() {
            assert_eq!(*record, expected_res[idx]);
        }
    }
}
