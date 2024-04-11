use clap::Parser;
use hic_lib::{read_data, write_rows_to_csv, ConverterError, HasChromosome, UniversalInteraction};
use serde::{Deserialize, Serialize};
use std::{
    cmp::max,
    cmp::min,
    collections::{HashMap, HashSet},
    error::Error,
    fs::{self, DirEntry, File, ReadDir},
    io::BufWriter,
    path::Path,
};

#[derive(Parser, Debug)]
struct Args {
    /// P value of files to read
    #[arg(short, long, default_value_t = 0.7)]
    pvalue: f32,

    /// Minimum intersection for two segments to be interpreted as intersecting
    #[arg(short, long, default_value_t = 0.2)]
    min_intersection: f64,

    #[arg(short, long, default_value = "../Processed_data/PCHi-C/Universal/")]
    interaction_path: String,

    #[arg(short, long)]
    hic_roadmap: Option<String>,

    #[arg(short, long)]
    gene_path: Option<String>,

    #[arg(short, long)]
    state_directory_path: Option<String>,

    #[arg(short, long, default_value = "../Processed_data/PCHi-C/Link_tissues/")]
    links_with_tissues: String,

    #[arg(
        short,
        long,
        default_value = "../Processed_data/PCHi-C/Interactions_with_states/"
    )]
    output_interaction_path: String,

    #[arg(short, long)]
    res_segment_lists: bool,
}

trait HasName {
    fn get_name(&self) -> &String;
}

trait HasStartEnd {
    fn get_start(&self) -> u32;
    fn get_end(&self) -> u32;
}

#[derive(Debug, Deserialize)]
struct GeneSegment {
    #[serde(alias = "name")]
    gene_name: String,
    #[serde(alias = "chrom")]
    chr: String,
    #[serde(alias = "txStart")]
    start: u32,
    #[serde(alias = "txEnd")]
    end: u32,
    #[serde(alias = "proteinID")]
    protein: String,
    #[serde(alias = "symbol")]
    gene_symbol: String,
}

impl HasName for GeneSegment {
    fn get_name(&self) -> &String {
        &self.gene_name
    }
}

impl HasStartEnd for GeneSegment {
    fn get_start(&self) -> u32 {
        self.start
    }

    fn get_end(&self) -> u32 {
        self.end
    }
}

#[derive(Debug, Deserialize)]
struct RoadmapHiCMapping {
    #[serde(alias = "Tissue_ID")]
    tissue_id: String,
    #[serde(alias = "Roadmap_ID")]
    roadmap_id: String,
}

#[derive(Debug, Deserialize)]
struct StateSegment {
    chr: String,
    start: u32,
    end: u32,
    state: String,
}

impl HasName for StateSegment {
    fn get_name(&self) -> &String {
        &self.state
    }
}

impl HasChromosome for StateSegment {
    fn get_chromosome(&self) -> &String {
        &self.chr
    }
}

impl HasStartEnd for StateSegment {
    fn get_start(&self) -> u32 {
        self.start
    }

    fn get_end(&self) -> u32 {
        self.end
    }
}

#[derive(Debug, Serialize)]
struct Link {
    bin1: u32,
    bin2: u32,
    bits: u32,
}

#[derive(Debug, Serialize)]
struct GeneSymbolAndProtein {
    symbol: String,
    protein: String,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct FinalData<'a> {
    tissue_data: &'a TissueData,
    states: Vec<String>,
    segments: Vec<SegmentWGenesStates>,
    links: Vec<Link>,
    gene_symbols_proteins: HashMap<String, GeneSymbolAndProtein>,
}

// Similar to FinalData but segments are not maps, instead there is one vector for segment start-end, segment genes, and segment states
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct FinalDataList<'a> {
    chr_names: Vec<String>,
    tissue_data: &'a TissueData,
    states: Vec<String>,
    chr_values: HashMap<String, ChrValues>,
    gene_symbols_proteins: HashMap<String, GeneSymbolAndProtein>,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct ChrValues {
    segments: Vec<Segment>,
    segment_genes: Vec<Vec<String>>,
    segment_states: Vec<HashMap<String, HashSet<String>>>,
    links: Vec<(u32, u32, u32)>,
}

#[derive(Debug, Serialize)]
struct SegmentWGenesStates {
    segment: Segment,
    genes: Vec<String>,
    states: HashMap<String, HashSet<String>>,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct TissueData {
    tissue_bits: HashMap<String, u32>,
    tissue_names: HashMap<String, String>,
    tissues_with_states: Vec<String>,
}

impl HasChromosome for GeneSegment {
    fn get_chromosome(&self) -> &String {
        &self.chr
    }
}

type Segment = (u32, u32);

fn get_state_to_category() -> HashMap<String, String> {
    let state_to_category: HashMap<String, String> = HashMap::from_iter(
        [
            ("1_TssA", "TSS"),
            ("2_TssFlnk", "TSS"),
            ("3_TssFlnkU", "TSS"),
            ("4_TssFlnkD", "TSS"),
            ("5_Tx", "TX"),
            ("6_TxWk", "TX"),
            ("7_EnhG1", "ENH"),
            ("8_EnhG2", "ENH"),
            ("9_EnhA1", "ENH"),
            ("10_EnhA2", "ENH"),
            ("11_EnhWk", "ENH"),
            ("12_ZNF/Rpts", "ZNF/RPTS"),
            ("13_Het", "HET"),
            ("14_TssBiv", "TSS_BIV"),
            ("15_EnhBiv", "ENH_BIV"),
            ("16_ReprPC", "PC"),
            ("17_ReprPCWk", "PC"),
            ("18_Quies", "QUI"),
        ]
        .iter()
        .map(|(k, v)| (k.to_string(), v.to_string())),
    );
    state_to_category
}

// Returns error if path does not exist
fn path_exists(input_path: &str) -> Result<(), Box<dyn Error>> {
    if !std::path::Path::new(&input_path).exists() {
        return Err(Box::new(ConverterError::FileError(format!(
            "input path '{}' does not exist",
            input_path
        ))));
    }
    Ok(())
}

// Gets paths of interaction data
fn get_interaction_data_paths(
    interaction_path: &String,
    pvalue: f32,
) -> Result<ReadDir, Box<dyn Error>> {
    let interaction_input_path = format!("{interaction_path}/Pvalue{pvalue}/");
    path_exists(&interaction_input_path)?;
    let paths = fs::read_dir(interaction_input_path)?;
    Ok(paths)
}

// Reads mapping of interaction data to state data
fn read_hic_to_roadmap(
    mapping_input_path: &String,
) -> Result<HashMap<String, String>, Box<dyn Error>> {
    path_exists(&mapping_input_path)?;
    let mut hic_roadmap: Vec<RoadmapHiCMapping> = Vec::new();
    read_data(&mut hic_roadmap, &mapping_input_path, b',', None)?;
    let hic_to_roadmap_mapping = HashMap::from_iter(
        hic_roadmap
            .iter()
            .map(|rec| (rec.tissue_id.clone(), rec.roadmap_id.clone())),
    );
    Ok(hic_to_roadmap_mapping)
}

// Gets sorted gene data, mapped to chromosomes
fn read_gene_data(
    genes_filename: &String,
) -> Result<HashMap<String, Vec<GeneSegment>>, Box<dyn Error>> {
    path_exists(&genes_filename)?;
    let mut gene_data: Vec<GeneSegment> = Vec::new();
    read_data(&mut gene_data, &genes_filename, b',', None)?;
    let mut chr_genes: HashMap<String, Vec<GeneSegment>> = map_to_chr(gene_data);
    // Sort genes by segment start
    for (_, value) in chr_genes.iter_mut() {
        value.sort_by_key(|k| k.start);
    }
    Ok(chr_genes)
}

// Creates output directory if it does not exist
fn create_output_directory(output_path: &String, pvalue: f32) -> Result<String, Box<dyn Error>> {
    let output_directory = format!("{output_path}Pvalue{pvalue}/")
        .to_string()
        .to_owned();
    fs::create_dir_all(&output_directory)?;
    Ok(output_directory)
}

fn path_to_filename(path: DirEntry) -> Result<String, Box<dyn Error>> {
    let filename = path.path();
    let filename = filename
        .to_str()
        .ok_or("Problem processing file path; check if path is valid unicode")?
        .to_string();
    Ok(filename)
}

// Reads interaction data and returns sorted segments and links, mapped to chromosomes
fn read_segments_links(
    filename: &String,
) -> Result<
    (
        HashMap<String, Vec<(u32, u32)>>,
        HashMap<String, Vec<((u32, u32), (u32, u32))>>,
    ),
    Box<dyn Error>,
> {
    // Read interaction data
    let mut interaction_data: Vec<UniversalInteraction> = Vec::new();
    read_data(&mut interaction_data, &filename, b',', None)?;
    let chr_interactions = map_to_chr(interaction_data);

    let chr_links: HashMap<String, Vec<(Segment, Segment)>> = get_links(&chr_interactions);
    let chr_segments: HashMap<String, Vec<Segment>> = get_segments(&chr_interactions);
    Ok((chr_segments, chr_links))
}

// Transforms a vector of data into a map, where each record is mapped to its corresponding chromosome
fn map_to_chr<T: HasChromosome>(data: Vec<T>) -> HashMap<String, Vec<T>> {
    let mut mapped_data: HashMap<String, Vec<T>> = HashMap::new();
    for record in data {
        let mut chr = record.get_chromosome().clone();
        if !chr.contains("chr") {
            chr.insert_str(0, "chr");
        }
        match mapped_data.get_mut(&chr) {
            Some(values) => {
                values.push(record);
            }
            None => {
                mapped_data.insert(chr.clone(), vec![record]);
            }
        };
    }
    mapped_data
}

// Gets a vector of sorted unique links from interactions for all chromocomes
fn get_links(
    chr_interactions: &HashMap<String, Vec<UniversalInteraction>>,
) -> HashMap<String, Vec<(Segment, Segment)>> {
    let mut chr_links: HashMap<String, Vec<(Segment, Segment)>> = HashMap::new();
    for (chr, interactions) in chr_interactions {
        // Create a set with all links in chromosome
        let mut links: HashSet<(Segment, Segment)> = HashSet::new();
        for interaction in interactions {
            let segm1: Segment = (interaction.bin1_start, interaction.bin1_end);
            let segm2: Segment = (interaction.bin2_start, interaction.bin2_end);
            links.insert((min(segm1, segm2), max(segm1, segm2)));
        }
        // Create a vector of segments, and sort it
        let mut links: Vec<(Segment, Segment)> = Vec::from_iter(links);
        links.sort();
        chr_links.insert(chr.clone(), links);
    }
    chr_links
}

// Gets a vector of unique sorted segments from interactions for all chromosomes
fn get_segments(
    chr_interactions: &HashMap<String, Vec<UniversalInteraction>>,
) -> HashMap<String, Vec<Segment>> {
    let mut chr_segments: HashMap<String, Vec<Segment>> = HashMap::new();
    for (chr, interactions) in chr_interactions {
        // Create a set with all segments in chromosome
        let mut segments: HashSet<Segment> = HashSet::new();
        for interaction in interactions {
            segments.extend(&[
                (interaction.bin1_start, interaction.bin1_end),
                (interaction.bin2_start, interaction.bin2_end),
            ]);
        }
        // Create a vector of segments, and sort it
        let mut segments = Vec::from_iter(segments);
        segments.sort();
        chr_segments.insert(chr.clone(), segments);
    }
    chr_segments
}

// Sets given bit for a list of links (chr_links) in all_links
fn set_link_tissue_bits(
    all_links: &mut HashMap<String, HashMap<(Segment, Segment), u32>>,
    chr_links: HashMap<String, Vec<(Segment, Segment)>>,
    bit: u32,
) {
    for chr in chr_links.keys() {
        if !all_links.contains_key(chr) {
            all_links.insert(chr.to_string(), HashMap::new());
        }
        let mutable_link_map = all_links.get_mut(chr).unwrap();
        // mutable_link_map.insert(((1, 2), (3, 4)), bit);
        for link in chr_links.get(chr).unwrap() {
            match mutable_link_map.get_mut(link) {
                Some(bits) => {
                    *bits |= bit;
                }
                None => {
                    mutable_link_map.insert(*link, bit);
                }
            }
        }
    }
}

// Updates set of all segments by adding new segments
fn add_new_segments(
    all_segments: &mut HashMap<String, HashSet<(u32, u32)>>,
    chr_segments: &HashMap<String, Vec<Segment>>,
) {
    for chr in chr_segments.keys() {
        match all_segments.get_mut(chr) {
            Some(segments) => {
                segments.extend(chr_segments.get(chr).unwrap().iter());
            }
            None => {
                let segments = chr_segments.get(chr).unwrap().clone();
                all_segments.insert(chr.clone(), HashSet::from_iter(segments));
            }
        }
    }
}

// Get tissue id from filename
fn get_tissue_id_and_name(filepath: &String) -> (String, String) {
    let filename = filepath.split("/").last().unwrap();
    let tissue_id = filename.split("_").next().unwrap();
    let prefix = tissue_id.to_string().clone() + "_";
    let tissue_name = filename
        .strip_prefix(prefix.as_str())
        .unwrap()
        .split(".")
        .next()
        .unwrap()
        .replace("_", " ");
    (tissue_id.to_string(), tissue_name.to_string())
}

// Returns a set of categories, into which states can be grouped
fn states_to_categories(states: Vec<String>) -> HashSet<String> {
    let state_categories: HashSet<String> = HashSet::from_iter(states.iter().map(|state| {
        get_state_to_category()
            .get(state)
            .expect(&format!("State with no category: {}", state))
            .clone()
    }));
    state_categories
}

fn read_state_data(
    state_directory_path: &String,
    roadmap_id: &String,
) -> Result<Vec<StateSegment>, Box<dyn Error>> {
    // Read state data
    let state_filename = state_directory_path.to_string() + "/" + roadmap_id + ".bed.gz";
    let mut state_data: Vec<StateSegment> = Vec::new();
    read_data(
        &mut state_data,
        &state_filename,
        b'\t',
        Some(vec!["chr", "start", "end", "state"]),
    )?;
    Ok(state_data)
}

// Reads states from specified Roadmap file, and adds states to segments in all_segment_states, if they intersect enough with state segments
fn update_segment_states(
    all_segment_states: &mut HashMap<String, HashMap<Segment, HashMap<String, HashSet<String>>>>,
    tissue_id: &String,
    chr_segments: &HashMap<String, Vec<(u32, u32)>>,
    min_intersection: f64,
    state_data: Vec<StateSegment>,
) {
    let chr_states: HashMap<String, Vec<StateSegment>> = map_to_chr(state_data);

    for chr in chr_segments.keys() {
        // Get a mapping of states to segments
        let segment_states: HashMap<Segment, Vec<String>> = map_properties_to_segments(
            chr_states.get(chr).unwrap(),
            chr_segments.get(chr).unwrap(),
            min_intersection,
        )
        .unwrap();
        // Update all_segment_states
        for (segment, states) in segment_states {
            update_chr_states(all_segment_states, chr, segment, tissue_id, states);
        }
    }
}

// Add states to all_segment_states[chr][segment][tissue_id]
// If any key does not exist yet, create it
fn update_chr_states(
    all_segment_states: &mut HashMap<String, HashMap<Segment, HashMap<String, HashSet<String>>>>,
    chr: &String,
    segment: Segment,
    tissue_id: &String,
    states: Vec<String>,
) {
    if !all_segment_states.contains_key(chr) {
        all_segment_states.insert(chr.clone(), HashMap::new());
    }

    let all_segment_states_chr: &mut HashMap<Segment, HashMap<String, HashSet<String>>> =
        all_segment_states.get_mut(chr).unwrap();

    if !all_segment_states_chr.contains_key(&segment) {
        all_segment_states_chr.insert(segment, HashMap::new());
    }

    let all_segment_states_chr_segm: &mut HashMap<String, HashSet<String>> =
        all_segment_states_chr.get_mut(&segment).unwrap();

    if !all_segment_states_chr_segm.contains_key(tissue_id) {
        all_segment_states_chr_segm.insert(tissue_id.to_string(), HashSet::new());
    }

    all_segment_states_chr_segm
        .get_mut(tissue_id)
        .unwrap()
        .extend(states_to_categories(states));
}

#[derive(Debug, Clone, PartialEq)]
pub enum DataError {
    NotSortedError(String),
}

impl std::fmt::Display for DataError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self {
            DataError::NotSortedError(msg) => write!(f, "{}", msg),
        }
    }
}

impl Error for DataError {}

// Maps properties (genes, states, or other) to interaction segment indices
fn map_properties_to_segments<T: HasStartEnd + HasName>(
    properties: &Vec<T>,
    interaction_segments: &Vec<Segment>,
    min_intersection: f64,
) -> Result<HashMap<Segment, Vec<String>>, Box<dyn Error>> {
    // If property vector is empty, return HashMap where empty vectors are mapped to segments
    if properties.len() == 0 {
        return Ok(HashMap::from_iter(
            interaction_segments.iter().map(|s| (*s, Vec::new())),
        ));
    }
    // Check if vectors are sorted
    if properties
        .windows(2)
        .all(|w| w[0].get_start() <= w[1].get_start())
        == false
    {
        return Err(Box::from(DataError::NotSortedError(
            "Property vector is not sorted by segment start".to_string(),
        )));
    }
    if interaction_segments.windows(2).all(|w| w[0].0 <= w[1].0) == false {
        return Err(Box::from(DataError::NotSortedError(
            "Interaction segment vector is not sorted by segment start".to_string(),
        )));
    }

    let mut relevant_properties: Vec<&T> = Vec::new();
    let mut segment_properties: HashMap<Segment, Vec<String>> = HashMap::new();
    let mut chr_properties_iter = properties.into_iter();
    let mut curr_property = chr_properties_iter.next().unwrap();

    // Iterate through interaction segments
    for segment in interaction_segments.iter() {
        segment_properties.insert(*segment, Vec::new());

        // Add properties whose segments potentially intersect with interaction segments to relevant_properties
        while curr_property.get_start() <= segment.1 {
            relevant_properties.push(&curr_property);
            match chr_properties_iter.next() {
                Some(property) => {
                    curr_property = property;
                }
                None => {
                    break;
                }
            }
        }

        // Create vector for properties that are still relevant
        let mut filtered_properties = Vec::new();
        for property in relevant_properties {
            // Add property to segment if segments intersect enough
            if intersects_enough(
                (property.get_start(), property.get_end()),
                (segment.0, segment.1),
                min_intersection,
            ) {
                segment_properties
                    .get_mut(segment)
                    .unwrap()
                    .push(property.get_name().clone());
            }
            // Keep properties that are still relevant
            if property.get_end() >= segment.0 {
                filtered_properties.push(property);
            }
        }
        relevant_properties = filtered_properties;
    }
    Ok(segment_properties)
}

// Checks if ratio (length of intersection of gene segment (gs) and interaction segment (is)) / (length of shortest segment of the two) is at least min_interaction
fn intersects_enough(ps: Segment, is: Segment, min_intersection: f64) -> bool {
    let intersection_len = min(ps.1, is.1).checked_sub(max(ps.0, is.0));
    match intersection_len {
        Some(len) => {
            return len as f64 / min(ps.1 - ps.0, is.1 - is.0) as f64 >= min_intersection;
            // return len as f64 / (ps.1 - ps.0) as f64 >= min_intersection;
        }
        None => {
            return false;
        }
    };
}

// Transforms links into sorted vector of Link structs where bin1 and bin2 are indices of segments
fn get_vec_of_links(
    links: &HashMap<(Segment, Segment), u32>,
    segment_to_idx: HashMap<Segment, usize>,
) -> Vec<Link> {
    let mut links: Vec<Link> = Vec::from_iter(links.iter().map(|link| Link {
        bin1: segment_to_idx[&link.0 .0].try_into().unwrap(),
        bin2: segment_to_idx[&link.0 .1].try_into().unwrap(),
        bits: *link.1,
    }));
    links.sort_unstable_by_key(|link| (link.bin1, link.bin2));
    links
}

// Create a file where each row represents a link (from id, to id, and a list of 0 and 1 (one value for each tissue), where 1 represents that link is present in tissue)
fn create_file_for_c(
    links: &Vec<Link>,
    tissues: &Vec<String>,
    chr: &String,
    links_with_tissues: &String,
    pvalue: f32,
) -> Result<(), Box<dyn Error>> {
    let mut links_for_c: Vec<Vec<u32>> = Vec::new();
    let mut header = vec!["chr".to_string(), "bait".to_string(), "oe".to_string()];
    header.append(&mut tissues.clone());
    for link in links.iter() {
        let chr_num = chr.strip_prefix("chr").unwrap_or(chr);
        let chr_num = chr_num.parse::<u32>().unwrap_or(23);
        let mut record = vec![chr_num, link.bin1, link.bin2];
        for (idx, _) in tissues.iter().enumerate() {
            if link.bits & 2_u32.pow(idx.try_into().unwrap()) > 0 {
                record.push(1);
            } else {
                record.push(0);
            }
        }
        links_for_c.push(record);
    }
    // Write data to a file
    write_rows_to_csv(
        &links_for_c,
        Some(header),
        format!("{links_with_tissues}/Pvalue{pvalue}/").as_str(),
        &(chr.clone() + ".csv"),
    )?;
    Ok(())
}

// Map symbols and proteins to gene names
fn map_genes_to_symbols_and_proteins(
    genes: &Vec<GeneSegment>,
) -> HashMap<String, GeneSymbolAndProtein> {
    let mut gene_symbols_proteins: HashMap<String, GeneSymbolAndProtein> = HashMap::new();
    for gene in genes {
        gene_symbols_proteins.insert(
            gene.gene_name.clone(),
            GeneSymbolAndProtein {
                symbol: gene.gene_symbol.clone(),
                protein: gene.protein.clone(),
            },
        );
    }
    gene_symbols_proteins
}

// Transforms data into preferred format for writing to file
fn get_final_results<'a>(
    segments: Vec<Segment>,
    all_links: &HashMap<(Segment, Segment), u32>,
    segment_genes: &mut HashMap<(u32, u32), Vec<String>>,
    gene_symbols_proteins: HashMap<String, GeneSymbolAndProtein>,
    all_segment_states: &mut HashMap<Segment, HashMap<String, HashSet<String>>>,
    tissue_data: &'a TissueData,
) -> FinalData<'a> {
    // Transforms links into vectors where each link has indices of two segments
    let segment_to_idx: HashMap<(u32, u32), usize> =
        HashMap::from_iter(segments.iter().enumerate().map(|(idx, &segm)| (segm, idx)));
    let links = get_vec_of_links(all_links, segment_to_idx);

    let segments_with_genes_states =
        Vec::from_iter(segments.iter().map(|&segm| SegmentWGenesStates {
            segment: segm,
            genes: segment_genes.get(&segm).unwrap().to_vec(),
            states: all_segment_states.remove(&segm).unwrap_or(HashMap::new()),
        }));

    let mut states: Vec<String> = Vec::from_iter(get_state_to_category().into_values());
    states.sort();
    states.dedup();

    let res: FinalData = FinalData {
        tissue_data: tissue_data.clone(),
        states,
        segments: segments_with_genes_states,
        links,
        gene_symbols_proteins,
    };
    res
}

// Creates FinalDataList from FinalData
fn final_data_to_list(final_data: FinalData, chr: String) -> FinalDataList {
    let mut segments: Vec<(u32, u32)> = Vec::new();
    let mut segment_genes: Vec<Vec<String>> = Vec::new();
    let mut segment_states: Vec<HashMap<String, HashSet<String>>> = Vec::new();
    for segm in final_data.segments {
        segments.push(segm.segment);
        segment_genes.push(segm.genes);
        segment_states.push(segm.states);
    }
    let links = Vec::from_iter(
        final_data
            .links
            .into_iter()
            .map(|x| (x.bin1, x.bin2, x.bits)),
    );
    let values: ChrValues = ChrValues {
        segments,
        segment_genes,
        segment_states,
        links,
    };
    let mut chr_values: HashMap<String, ChrValues> = HashMap::new();
    chr_values.insert(chr.clone(), values);
    let final_data_list: FinalDataList = FinalDataList {
        chr_names: vec![chr],
        tissue_data: final_data.tissue_data,
        states: final_data.states,
        chr_values,
        gene_symbols_proteins: final_data.gene_symbols_proteins,
    };
    final_data_list
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    // Read everything and prepare output directory
    let paths = get_interaction_data_paths(&args.interaction_path, args.pvalue)?;
    let states_available = if let (Some(hic_roadmap), Some(state_directory)) = (
        args.hic_roadmap.as_ref(),
        args.state_directory_path.as_ref(),
    ) {
        Path::new(hic_roadmap).exists() && Path::new(state_directory).exists()
    } else {
        false
    };
    if !states_available {
        println!("Mapping of state to interaction tissue names or state directory not available")
    }
    let mut chr_genes = HashMap::new();
    let mut genes_available = false;
    if let Some(gene_path) = args.gene_path.as_ref() {
        genes_available = Path::new(gene_path).exists();
        if genes_available {
            chr_genes = read_gene_data(gene_path)?;
        }
    }
    if !genes_available {
        println!("Genes not available")
    }
    let output_directory = create_output_directory(&args.output_interaction_path, args.pvalue)?;

    // Create empty containers for data across all tissues
    let mut all_segments: HashMap<String, HashSet<Segment>> = HashMap::new();
    let mut all_links: HashMap<String, HashMap<(Segment, Segment), u32>> = HashMap::new();
    let mut all_segment_states: HashMap<
        String,
        HashMap<Segment, HashMap<String, HashSet<String>>>,
    > = HashMap::new();
    let mut tissue_data: TissueData = TissueData {
        tissue_bits: HashMap::new(),
        tissue_names: HashMap::new(),
        tissues_with_states: Vec::new(),
    };
    let mut tissues: Vec<String> = Vec::new();

    // Iterate through interaction data files
    for (tissue_idx, path) in paths.enumerate() {
        let filename = path_to_filename(path?)?;

        let (chr_segments, chr_links) = read_segments_links(&filename)?;

        // Update all_links and all_segments with data of current tissue
        let bit = 2_u32.pow(tissue_idx.try_into().unwrap());
        set_link_tissue_bits(&mut all_links, chr_links, bit);
        add_new_segments(&mut all_segments, &chr_segments);

        let (tissue_id, tissue_name) = get_tissue_id_and_name(&filename);
        tissue_data
            .tissue_names
            .insert(tissue_id.clone(), tissue_name);
        tissues.push(tissue_id.clone());
        tissue_data.tissue_bits.insert(tissue_id.clone(), bit);

        if let (Some(hic_roadmap), Some(state_directory)) = (
            args.hic_roadmap.as_ref(),
            args.state_directory_path.as_ref(),
        ) {
            if states_available {
                // Update states if they are available for current tissue
                let hic_to_roadmap_mapping = read_hic_to_roadmap(hic_roadmap)?;
                let roadmap_id = hic_to_roadmap_mapping.get(&tissue_id);

                if path_exists(state_directory).is_ok()
                    && roadmap_id.is_some()
                    && roadmap_id.unwrap() != ""
                {
                    let state_data = read_state_data(state_directory, roadmap_id.unwrap());
                    // let x = state_data.as_ref().unwrap();
                    if state_data.is_ok() {
                        update_segment_states(
                            &mut all_segment_states,
                            &tissue_id,
                            &chr_segments,
                            args.min_intersection,
                            state_data.unwrap(),
                        );
                    }
                    tissue_data.tissues_with_states.push(tissue_id.clone());
                }
            }
        }
    }

    for chr in all_segments.keys() {
        // Get sorted lists of all segments and links in current chromosome
        let mut segments: Vec<Segment> = Vec::from_iter(all_segments.get(chr).unwrap().clone());
        segments.sort();

        let mut gene_symbols_proteins = HashMap::new();
        let mut genes = &Vec::new();
        if genes_available {
            genes = chr_genes.get(chr).unwrap();
            // Map gene names to symbols and proteins
            gene_symbols_proteins = map_genes_to_symbols_and_proteins(genes);
        }
        // Map genes to interaction segment indices
        let mut segment_genes =
            map_properties_to_segments(genes, &segments, args.min_intersection).unwrap();

        let res: FinalData = get_final_results(
            segments,
            all_links.get(chr).unwrap(),
            &mut segment_genes,
            gene_symbols_proteins,
            all_segment_states
                .get_mut(chr)
                .unwrap_or(&mut HashMap::new()),
            &tissue_data,
        );

        // Write links and tissues to csv file for component finder program
        create_file_for_c(
            &res.links,
            &tissues,
            chr,
            &args.links_with_tissues,
            args.pvalue,
        )?;

        // Write final data to file
        if !std::path::Path::new(&output_directory).exists() {
            return Err(Box::new(ConverterError::FileError(format!(
                "output directory '{}' does not exist",
                output_directory
            ))));
        }
        let output_path = output_directory.to_string() + chr + ".json";
        let output_file = File::create(output_path).expect("Could not open output file");
        let writer = BufWriter::new(output_file);

        if args.res_segment_lists {
            let res = final_data_to_list(res, chr.to_string());
            serde_json::to_writer(writer, &res).expect("Failed serializing and writing results");
        } else {
            serde_json::to_writer(writer, &res).expect("Failed serializing and writing results");
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_intersects_enough() {
        let e7 = 10_u32.pow(7);
        // Barely intersecting enough
        assert_eq!(
            intersects_enough((1 * e7, 5 * e7), (3 * e7, 13 * e7), 0.2),
            true
        );
        // Barely not intersecting enough
        assert_eq!(
            intersects_enough((1 * e7, 5 * e7), (3 * e7 + 1, 13 * e7), 0.2),
            false
        );
        // One segment is entirely inside the other
        assert_eq!(
            intersects_enough((5 * e7, 6 * e7), (4 * e7, 7 * e7), 0.2),
            true
        );
        // One small segment entirely inside the other
        assert_eq!(
            intersects_enough((3 * e7, 4 * e7), (1 * e7, 10 * e7), 0.2),
            false
        );
        // Segments one after another, property segment first
        assert_eq!(
            intersects_enough((5 * e7, 6 * e7), (6 * e7, 8 * e7), 0.2),
            false
        );
        // Segments one after another, interaction segment first
        assert_eq!(
            intersects_enough((5 * e7, 6 * e7), (6 * e7, 8 * e7), 0.2),
            false
        );
        // Barely intersecting enough
        assert_eq!(
            intersects_enough((1 * e7, 5 * e7), (3 * e7 + 100, 13 * e7), 0.1),
            true
        );
    }

    fn assert_my_property_maps_equal(
        res: HashMap<Segment, Vec<String>>,
        correct_mapping: HashMap<Segment, Vec<String>>,
    ) {
        // Check if segment sets are equal
        let res_segments: HashSet<&(u32, u32)> = HashSet::from_iter(res.keys().into_iter());
        let expected_segments: HashSet<&(u32, u32)> =
            HashSet::from_iter(correct_mapping.keys().into_iter());
        assert_eq!(res_segments, expected_segments);

        // Check if states for each segment are equal
        for segm in expected_segments {
            let mut res_states = res.get(segm).unwrap().clone();
            res_states.sort();
            let mut expected_states = correct_mapping.get(segm).unwrap().clone();
            expected_states.sort();
            assert_eq!(res_states, expected_states);
        }
    }
    #[test]
    fn test_map_properties_to_segments() {
        let min_intersection = 0.2;
        let mut properties: Vec<StateSegment> = Vec::new();
        let mut interaction_segments: Vec<Segment> = Vec::new();
        let mut correct_mapping: HashMap<Segment, Vec<String>> = HashMap::new();

        for i in 0..20 {
            let segm = (2 * i * 117, (2 * i + 1) * 117);
            interaction_segments.push(segm);
            correct_mapping.insert(segm, Vec::new());
        }
        // Test with empty vector of states
        let res = map_properties_to_segments(&properties, &interaction_segments, min_intersection)
            .unwrap();
        assert_eq!(res, correct_mapping);

        // Function does not use chr so it can be any string
        let chr = "chr13".to_string();
        // Normal test with one state
        let state = "18_Qui".to_string();
        let property_segments: Vec<Segment> = vec![
            (0, 171),
            (494, 513),
            (684, 741),
            (1026, 1197),
            (1482, 1539),
            (1710, 1729),
            (2052, 2223),
            (2470, 2565),
            (3078, 3211),
            (3458, 3591),
            (4104, 4199),
        ];
        for segm in property_segments {
            properties.push(StateSegment {
                chr: chr.clone(),
                start: segm.0,
                end: segm.1,
                state: state.clone(),
            });
        }
        let segments_with_qui = vec![
            (0, 117),
            (702, 819),
            (936, 1053),
            (1170, 1287),
            (1404, 1521),
            (2106, 2223),
            (3042, 3159),
            (3510, 3627),
        ];
        for segm in segments_with_qui {
            correct_mapping.get_mut(&segm).unwrap().push(state.clone());
        }
        let res = map_properties_to_segments(&properties, &interaction_segments, min_intersection)
            .unwrap();
        assert_eq!(res, correct_mapping);

        // Add one different state to end
        properties.push(StateSegment {
            chr: chr.clone(),
            start: 4446,
            end: 4617,
            state: "6_TxWk".to_string(),
        });
        correct_mapping
            .get_mut(&(4446, 4563))
            .unwrap()
            .push("6_TxWk".to_string());
        let res = map_properties_to_segments(&properties, &interaction_segments, min_intersection)
            .unwrap();
        assert_eq!(res, correct_mapping);

        // Add one different state somewhere in the middle without sorting
        properties.push(StateSegment {
            chr: chr.clone(),
            start: 900,
            end: 1000,
            state: "5_Tx".to_string(),
        });
        correct_mapping
            .get_mut(&(936, 1053))
            .unwrap()
            .push("5_Tx".to_string());
        let res = map_properties_to_segments(&properties, &interaction_segments, min_intersection);
        // Since segments are no longer sorted, expect an error
        assert!(res.is_err());

        properties.sort_unstable_by_key(|prop| prop.get_start());
        let res = map_properties_to_segments(&properties, &interaction_segments, min_intersection)
            .unwrap();
        // Since order of states in vectors does not matter, use custom equality assertion
        assert_my_property_maps_equal(res, correct_mapping);
    }
}
