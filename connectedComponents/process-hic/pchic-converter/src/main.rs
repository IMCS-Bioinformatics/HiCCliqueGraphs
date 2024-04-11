use clap::Parser;
use hic_lib::{read_data, write_structs_to_csv, ConverterError, UniversalInteraction};
use serde::Deserialize;
use std::{collections::HashMap, error::Error, fs};

#[derive(Parser, Debug)]
struct Args {
    /// P value by which to filter data
    #[arg(short, long, default_value = "0.7")]
    pvalue: f32,

    #[arg(short, long, default_value = "../Source_data/pcHi-C(hg19)/")]
    og_interaction_path: String,

    #[arg(
        short,
        long,
        default_value = "../Source_data/pcHi-C(hg19)/tissue-id-name-roadmap.csv"
    )]
    tissue_id_names_path: String,

    #[arg(short, long, default_value = "../Processed_data/PCHiC/Universal/")]
    universal_data_path: String,
}

#[derive(Debug, Deserialize)]
struct TissueNameID {
    #[serde(alias = "Tissue_name")]
    tissue_name: String,
    #[serde(alias = "Tissue_ID")]
    tissue_id: String,
}

// Interaction contains the main fields from original PCHi-C interaction data.
#[derive(Debug, Deserialize)]
struct Interaction {
    frag1: String,
    frag2: String,
    #[serde(alias = "-log10(result)")]
    pvalue: f32,
}

// Returns interaction transformed into universal format
fn make_universal_record(rec: &Interaction) -> UniversalInteraction {
    let bin1 = rec.frag1.split(':').last().unwrap();
    let bin2 = rec.frag2.split(':').last().unwrap();
    let universal_record: UniversalInteraction = UniversalInteraction {
        chr: rec.frag1.split(':').nth(0).unwrap().to_string(),
        bin1_start: bin1.split('-').nth(0).unwrap().parse().unwrap(),
        bin1_end: bin1.split('-').nth(1).unwrap().parse().unwrap(),
        bin2_start: bin2.split('-').nth(0).unwrap().parse().unwrap(),
        bin2_end: bin2.split('-').nth(1).unwrap().parse().unwrap(),
        p_value: rec.pvalue,
    };
    universal_record
}

// Creates a sorted vector of universal format interaction data, filtering out non-significant interactions
// Interaction between same pairs of chromatin segments may still remain, and individual interactions are not sorted (bin1 may be greater than bin2)
fn transform_data_to_universal(data: &Vec<Interaction>, pvalue: f32) -> Vec<UniversalInteraction> {
    let mut universal_data = Vec::new();
    for record in data {
        if record.pvalue >= pvalue {
            if record.frag1.split(':').nth(0) == record.frag2.split(':').nth(0) {
                universal_data.push(make_universal_record(&record));
            } else {
                println!("{:?}", record);
            }
        }
    }
    universal_data.sort_by_key(|item| (item.chr.clone(), item.bin1_start, item.bin2_start));
    universal_data
}

fn create_output_filename(
    filename: &str,
    tissue_id_names: &HashMap<String, String>,
) -> Result<String, Box<dyn Error>> {
    let tissue_id: &str = filename
        .split('/')
        .last()
        .ok_or("could not get filename ending")?
        .split(".")
        .nth(0)
        .ok_or("could not get tissue id from filename")?;
    let output_filename = format!(
        "{}_{}.csv",
        tissue_id,
        tissue_id_names
            .get(tissue_id)
            .unwrap_or(&tissue_id.to_string())
    );
    Ok(output_filename)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    if !std::path::Path::new(&args.og_interaction_path).exists() {
        return Err(Box::new(ConverterError::FileError(format!(
            "input directory '{}' does not exist",
            &args.og_interaction_path
        ))));
    }
    let paths = fs::read_dir(&args.og_interaction_path)?;
    let output_directory = format!("{0}/Pvalue{1}/", args.universal_data_path, args.pvalue)
        .to_string()
        .to_owned();
    fs::create_dir_all(&output_directory)?;
    let mut tissue_id_names_vec: Vec<TissueNameID> = Vec::new();
    let mut tissue_id_names: HashMap<String, String> = HashMap::new();
    if read_data(
        &mut tissue_id_names_vec,
        &args.tissue_id_names_path,
        b',',
        None,
    )
    .is_ok()
    {
        tissue_id_names = HashMap::from_iter(
            tissue_id_names_vec
                .iter()
                .map(|rec| (rec.tissue_id.clone(), rec.tissue_name.clone())),
        );
    }
    for path in paths {
        let file_path = path?.path();
        let mut og_data: Vec<Interaction> = Vec::new();
        let mut og_po_data: Vec<Interaction> = Vec::new();
        let filename = file_path.to_str().ok_or(format!(
            "Problem processing file path; check if path {} is valid unicode",
            file_path.display()
        ))?;
        if !filename.ends_with("pp.txt.zip") {
            continue;
        }
        read_data(&mut og_data, &filename, b'\t', None)?;
        read_data(&mut og_po_data, &filename.replace("pp", "po"), b'\t', None)?;
        og_data.append(&mut og_po_data);
        let universal_data = transform_data_to_universal(&og_data, args.pvalue);
        write_structs_to_csv(
            &universal_data,
            &output_directory,
            &create_output_filename(filename, &tissue_id_names)?,
        )?;
    }
    Ok(())
}
