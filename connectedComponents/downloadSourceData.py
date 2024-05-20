import urllib.request
import os

replace_existing_files=False

def download_file(url, local_path):
    if not replace_existing_files and os.path.exists(local_path):
        raise Exception(f"Local file already exists: {os.getcwd()}/{local_path}")
    else:
        print(f"Downloading file: {url}", flush=True)
        urllib.request.urlretrieve(url, local_path)

os.chdir(os.path.dirname(__file__))

data_url = "http://susurs.mii.lu.lv/HiCData/"
roadmap_ids = ["E080", "E065", "E066", "E096", "E097", "E098", "E100", "E109", "E113"]
hic_tissue_names_ids = {"Adrenal_gland": "AD", "Aorta": "AO", "Bladder": "BL", "Liver": "LI", "Lung": "LG", "Ovary": "OV", "Pancreas": "PA", "Psoas": "PO", "Small_Bowel": "SB", "Spleen": "SX"}
pchic_tissue_ids = ["AD2", "AO", "BL1", "LG", "LI11", "OV2", "PA", "PO3", "SB", "SX"]

state_path = "Source_data/chromatin_states/"
state_path_bed = f"{state_path}18states.all.mnemonics.bedFiles/"
os.makedirs(state_path_bed, exist_ok=True)
normal_hic_path = "Source_data/Normal_HiC(hg19,3DIV_legacy)_10/"
os.makedirs(normal_hic_path, exist_ok=True)
pchic_path = "Source_data/pcHi-C(hg19)/"
os.makedirs(pchic_path, exist_ok=True)
tad_path = "Processed_data/NormalHiC/TADs/"
os.makedirs(tad_path, exist_ok=True)

for id in roadmap_ids:
    download_file(f"{data_url}18states.all.mnemonics.bedFiles/{id}_18_core_K27ac_mnemonics.bed.gz", f"{state_path_bed}{id}.bed.gz")

download_file(f"{data_url}Tissue_name_mapping/roadmap-normalHiC.csv", f"{state_path}roadmap-normalHiC.csv")
download_file(f"{data_url}Tissue_name_mapping/roadmap-pcHiC.csv", f"{state_path}roadmap-pcHiC.csv")

for name, id in hic_tissue_names_ids.items():
    os.makedirs(f"{normal_hic_path}{name}", exist_ok=True)
    download_file(f"{data_url}Normal_HiC(hg19,3DIV_legacy)_selection/{name}/3DIV_{id}_cutoff-10.zip", f"{normal_hic_path}{name}/3DIV_{id}_cutoff-10.zip")
    download_file(f"{data_url}Normal_HiC(hg19,3DIV_legacy)_selection/{name}/3DIV_{id}_cutoff-2.zip", f"{normal_hic_path}{name}/3DIV_{id}_cutoff-2.zip")

for id in pchic_tissue_ids:
    download_file(f"{data_url}pcHi-C(hg19)/{id}.po.txt.zip", f"{pchic_path}{id}.po.txt.zip")
    download_file(f"{data_url}pcHi-C(hg19)/{id}.pp.txt.zip", f"{pchic_path}{id}.pp.txt.zip")

download_file(f"{data_url}Tissue_name_mapping/tissue-id-name-roadmap.csv", f"{pchic_path}tissue-id-name-roadmap.csv")

download_file(f"{data_url}hg19_genes_symbols.csv", "Source_data/hg19_genes+symbols.csv")

chrs = [i for i in range(1, 23)] + ["X"]
for chr in chrs:
    download_file(f"{data_url}TADs/RealRes-chr{chr}-SpectralTad.json", f"{tad_path}RealRes-chr{chr}-SpectralTad.json")
