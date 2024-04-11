#!/bin/bash
dataPath="Processed_data/PCHiC/"
pvalue=0.7
trap "exit" INT
# Run Rust
cd process-hic
cargo run --bin pchic-converter --release
cargo run --bin add-states-and-genes --release -- --pvalue $pvalue --interaction-path "../${dataPath}Universal" --hic-roadmap "../Source_data/chromatin_states/roadmap-pcHiC.csv" --gene-path "../Source_data/hg19_genes+symbols.csv" --state-directory-path "../Source_data/chromatin_states/18states.all.mnemonics.bedFiles" --links-with-tissues "../${dataPath}Link_tissues/" --output-interaction-path "../${dataPath}Interactions_with_states/"
cd ../

# Run component finder
g++ -m32 -O3 -std=c++11 Find_Components/FindNetworkComponents.cpp -o Find_Components/FindNetworkComponents.exe
inputPath="${dataPath}Link_tissues/Pvalue${pvalue}"
outputPath="${dataPath}Components/Pvalue${pvalue}"

mkdir -p ${outputPath}
for n in {1..22}
do
	"Find_Components/FindNetworkComponents.exe" "${inputPath}/chr"$n".csv" "${outputPath}/chr"$n"_result.txt" "${outputPath}/chr"$n".csv" "100" "10" "0.75" "0.25"
done
"Find_Components/FindNetworkComponents.exe" "${inputPath}/chrX.csv" "${outputPath}/chrX_result.txt" "${outputPath}/chrX.csv" "100" "10" "0.75" "0.25"

#Run Python programs
python Python_programs/calcStats.py --dataPath $dataPath --pvalue $pvalue
python Python_programs/getStatePercentageFull.py --dataPath $dataPath --pvalue $pvalue
python Python_programs/makeZipFiles.py --dataPath $dataPath --pvalue $pvalue