use csv::{ReaderBuilder, StringRecord};
use flate2::read::GzDecoder;
use serde::{de, Deserialize, Serialize};
use std::{
    error::Error,
    fs,
    io::{self},
};

// UniversalInteraction contains slightly edited Hi-C interaction data in universal format.
// Chromatin segments have a start and an end.
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct UniversalInteraction {
    pub chr: String,
    pub bin1_start: u32,
    pub bin1_end: u32,
    pub bin2_start: u32,
    pub bin2_end: u32,
    pub p_value: f32,
}

pub trait HasChromosome {
    fn get_chromosome(&self) -> &String;
}

impl HasChromosome for UniversalInteraction {
    fn get_chromosome(&self) -> &String {
        &self.chr
    }
}

#[derive(Debug, Clone)]
pub enum ConverterError {
    FileError(String),
}

impl std::fmt::Display for ConverterError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self {
            ConverterError::FileError(msg) => write!(f, "{}", msg),
        }
    }
}

impl Error for ConverterError {}

// Reads data from path filename into data; file may be compressed (.zip)
pub fn read_data<T: de::DeserializeOwned>(
    data: &mut Vec<T>,
    filename: &str,
    delimiter: u8,
    header: Option<Vec<&str>>,
) -> Result<(), Box<dyn Error>> {
    let file_metadata = fs::metadata(&filename)?;
    if !file_metadata.is_file() {
        return Err(Box::from(ConverterError::FileError(format!(
            "'{}' is not a file",
            filename
        ))));
    }
    let read_file = || -> Result<(), Box<dyn Error>> {
        let f = std::fs::File::open(&filename)?;
        if filename.ends_with("zip") {
            let mut ar = zip::read::ZipArchive::new(&f)?;
            let zf = ar.by_index(0)?;
            read_dsv(data, delimiter, zf, header)?;
        } else if filename.ends_with("gz") {
            let gz = GzDecoder::new(f);
            read_dsv(data, delimiter, gz, header)?;
        } else {
            read_dsv(data, delimiter, f, header)?;
        };
        Ok(())
    };
    if let Err(err) = read_file() {
        return Err(Box::from(ConverterError::FileError(format!(
            "Problem reading file '{}', {}",
            filename, err
        ))));
    }
    Ok(())
}

fn read_dsv<TData: de::DeserializeOwned, TReader: io::Read>(
    data: &mut Vec<TData>,
    delimiter: u8,
    reader: TReader,
    header: Option<Vec<&str>>,
) -> Result<(), Box<dyn Error>> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(delimiter)
        .from_reader(reader);
    if header.is_some() {
        rdr.set_headers(StringRecord::from(header.unwrap()));
    }
    for result in rdr.deserialize() {
        let record: TData = result?;
        data.push(record);
    }
    Ok(())
}

// Writes serializable data to <filename> in <directory> as csv file
pub fn write_structs_to_csv<T: serde::Serialize>(
    data: &Vec<T>,
    directory: &str,
    filename: &str,
) -> Result<(), Box<dyn Error>> {
    fs::create_dir_all(&directory)?;
    let mut wtr = csv::Writer::from_path(directory.to_string() + filename)?;
    for record in data {
        wtr.serialize(record)?;
    }
    wtr.flush()?;
    Ok(())
}

// Writes rows of data to <filename> in <directory> as csv file
pub fn write_rows_to_csv<T: serde::Serialize>(
    data_rows: &Vec<Vec<T>>,
    header: Option<Vec<String>>,
    directory: &str,
    filename: &str,
) -> Result<(), Box<dyn Error>> {
    fs::create_dir_all(&directory)?;
    let mut wtr = csv::Writer::from_path(directory.to_string() + filename)?;
    if header.is_some() {
        wtr.write_record(header.unwrap())?;
    }
    for row in data_rows {
        wtr.serialize(row)?;
    }
    wtr.flush()?;
    Ok(())
}
