#!/bin/bash

# ============================================================================
# run_full_imputation.sh - Full dataset imputation pipeline
# ============================================================================
# Usage: ./run_full_imputation.sh [scool_file] [chromosomes] [r_values]
# Example: ./run_full_imputation.sh sourceData/nagano_10kb_cell_types.scool "chr1 chr2" "1.3 2.0"
# Example: ./run_full_imputation.sh  (uses defaults)
# ============================================================================

set -e  # Exit on error

# Configuration
SCOOL_FILE="${1:-sourceData/nagano_10kb_cell_types.scool}"
CHROMOSOMES="${2:-ALL}"  # Space-separated list or "ALL"
R_VALUES="${3:-1.3}"     # Space-separated list of r values (mean node degree)

# Imputation parameters
PAD=1
STD=1
RP=0.5

# Directories (script runs from imputation directory)
BASE_DIR="."
RAW_DIR="raw"
IMPUTED_DIR="imputed"
CHROM_SIZES_FILE="chrom_sizes.txt"

# Parallel processing
NUM_PROCESSES=4

echo ""
echo "=========================================="
echo "SCHiCRank Full Imputation Pipeline"
echo "=========================================="
echo "Input scool:  $SCOOL_FILE"
echo "Chromosomes:  $CHROMOSOMES"
echo ""
echo "Imputation parameters:"
echo "  PAD:          $PAD"
echo "  STD:          $STD"
echo "  RP:           $RP"
echo ""
echo "Sparsity levels (mean node degree):"
echo "  r values:     $R_VALUES"
echo ""
echo "=========================================="
echo ""

# Activate schicluster environment
echo "Activating schicluster environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
CONDA_ENV_NAME=$(grep '^CONDA_ENV_NAME=' ../setup_schicluster_env.sh | cut -d'"' -f2)
conda activate "$CONDA_ENV_NAME"

# Detect resolution from scool file
echo "Detecting resolution from scool file..."
RESOLUTION=$(python -c "import h5py, cooler; f = h5py.File('$SCOOL_FILE', 'r'); cell = list(f['cells'].keys())[0]; f.close(); print(cooler.Cooler('$SCOOL_FILE::cells/' + cell).binsize)")
echo "Resolution: ${RESOLUTION}bp"
echo ""

# Step 1: Generate chromosome sizes if not exists
if [ ! -f "$CHROM_SIZES_FILE" ]; then
    echo "📋 Generating chromosome sizes file..."
    python exportScoolToText.py \
        --scool-file "$SCOOL_FILE" \
        --chrom "chr1" \
        --output-dir "$RAW_DIR/chr1" \
        --chrom-sizes-file "$CHROM_SIZES_FILE" \
        --max-cells 1
    echo ""
fi

# Step 2: Get list of chromosomes to process
if [ "$CHROMOSOMES" == "ALL" ]; then
    CHROM_LIST=$(awk '{print $1}' "$CHROM_SIZES_FILE" | grep -E '^chr[0-9]+$|^chrX$')
else
    CHROM_LIST="$CHROMOSOMES"
fi

echo "📂 Processing chromosomes: $CHROM_LIST"
echo ""

# Step 3: Export cells from scool to text format
echo "=========================================="
echo "Step 1: Exporting cells from scool..."
echo "=========================================="

for CHROM in $CHROM_LIST; do
    # Check if already exported
    if [ -d "$RAW_DIR/$CHROM" ] && [ -n "$(ls -A $RAW_DIR/$CHROM/*.txt 2>/dev/null)" ]; then
        NUM_FILES=$(ls -1 $RAW_DIR/$CHROM/*.txt 2>/dev/null | wc -l)
        echo "✓ Chromosome $CHROM already exported ($NUM_FILES cells found)"
        echo "  To re-export, delete: $RAW_DIR/$CHROM"
        continue
    fi

    echo "Exporting chromosome: $CHROM"

    python exportScoolToText.py \
        --scool-file "$SCOOL_FILE" \
        --chrom "$CHROM" \
        --output-dir "$RAW_DIR/$CHROM"

    if [ $? -ne 0 ]; then
        echo "❌ Failed to export cells for $CHROM"
        exit 1
    fi
done

echo ""
echo "=========================================="
echo "Step 2: Running scHiCluster imputation..."
echo "=========================================="

# Step 4: Run imputation for each chromosome
for CHROM in $CHROM_LIST; do
    echo "🔄 Imputing chromosome: $CHROM"

    mkdir -p "$IMPUTED_DIR/$CHROM"

    # Get list of cell files
    CELL_FILES=("$RAW_DIR/$CHROM/"*_${CHROM}.txt)
    # Check if glob matched anything
    if [ ! -f "${CELL_FILES[0]}" ]; then
        echo "  No cells found for $CHROM, skipping."
        continue
    fi
    TOTAL_CELLS=${#CELL_FILES[@]}

    echo "Found $TOTAL_CELLS cells to impute for $CHROM"

    # Create chromosome-specific sizes file
    CHROM_SPECIFIC_FILE="${CHROM}_sizes.txt"
    if [ ! -f "$CHROM_SPECIFIC_FILE" ]; then
        grep "^${CHROM}\s" "$CHROM_SIZES_FILE" > "$CHROM_SPECIFIC_FILE"
    fi

    # Process cells
    CURRENT=0
    for CELL_FILE in "${CELL_FILES[@]}"; do
        CURRENT=$((CURRENT + 1))

        # Extract cell name from filename
        CELL_NAME=$(basename "$CELL_FILE" "_${CHROM}.txt")

        # Check if already imputed (look for output with std as float)
        EXISTING_OUTPUT=$(ls "$IMPUTED_DIR/$CHROM/${CELL_NAME}_${CHROM}_pad${PAD}_std"*"_rp${RP}_sqrtvc.hdf5" 2>/dev/null | head -1)
        if [ -f "$EXISTING_OUTPUT" ]; then
            if [ $((CURRENT % 50)) -eq 0 ]; then
                echo "  Cell $CURRENT/$TOTAL_CELLS: $CELL_NAME (already imputed)"
            fi
            continue
        fi

        if [ $((CURRENT % 10)) -eq 0 ]; then
            echo "  Processing cell $CURRENT/$TOTAL_CELLS: $CELL_NAME"
        fi

        # Run hicluster imputation
        hicluster impute-cell \
            --indir "$RAW_DIR/$CHROM/" \
            --outdir "$IMPUTED_DIR/$CHROM/" \
            --cell "$CELL_NAME" \
            --chrom "$CHROM" \
            --res "$RESOLUTION" \
            --chrom_file "$CHROM_SPECIFIC_FILE" \
            --pad "$PAD" \
            --std "$STD" \
            --rp "$RP" \
            > /dev/null 2>&1 || echo "  ⚠️  Warning: Imputation failed for $CELL_NAME"
    done

    echo "✅ Completed imputation for $CHROM"
done

echo ""
echo "=========================================="
echo "Step 3: Filtering by sparsity levels..."
echo "=========================================="

# Calculate total bins
CHROM_SIZES_TOTAL=0
for CHROM in $CHROM_LIST; do
    SIZE=$(awk -v chr="$CHROM" '$1==chr {print $2}' "$CHROM_SIZES_FILE")
    N_BINS=$(( (SIZE + RESOLUTION - 1) / RESOLUTION ))
    CHROM_SIZES_TOTAL=$(( CHROM_SIZES_TOTAL + N_BINS ))
done

echo "Total bins across chromosomes: $CHROM_SIZES_TOTAL"

# For each r value, filter the imputed data
for R in $R_VALUES; do
    echo ""
    echo "Processing r = $R (mean node degree)..."

    FILTERED_DIR="filtered_r${R}"
    mkdir -p "$FILTERED_DIR"

    # Filter each chromosome's imputed HDF5 files (batch mode - single Python process per chromosome)
    for CHROM in $CHROM_LIST; do
        echo "  Filtering chromosome: $CHROM"

        CHROM_SIZE=$(awk -v chr="$CHROM" '$1==chr {print $2}' "$CHROM_SIZES_FILE")
        N_BINS=$(( (CHROM_SIZE + RESOLUTION - 1) / RESOLUTION ))

        TOTAL_CELLS=$(ls "$IMPUTED_DIR/$CHROM"/*.hdf5 2>/dev/null | wc -l)
        if [ "$TOTAL_CELLS" -eq 0 ]; then
            echo "  No imputed files for $CHROM, skipping."
            continue
        fi
        echo "    $TOTAL_CELLS cells to filter"

        python filterImputedBySparsity.py \
            --input-dir "$IMPUTED_DIR/$CHROM" \
            --output-dir "$FILTERED_DIR/$CHROM" \
            --r "$R" \
            --n-bins "$N_BINS"

        if [ $? -ne 0 ]; then
            echo "    ⚠️  Warning: Filtering failed for $CHROM"
        fi
    done

    echo "  ✅ Filtering complete for r = $R"
done

echo ""
echo "=========================================="
echo "Step 4: Converting to scool format..."
echo "=========================================="

# Convert filtered HDF5 files to scool format for each r value
OUTPUT_SCOOLS=()

for R in $R_VALUES; do
    echo ""
    echo "Creating scool for r = $R..."

    FILTERED_DIR="filtered_r${R}"

    # Create filename with all imputation parameters
    SCOOL_BASENAME=$(basename "$SCOOL_FILE" .scool)
    OUTPUT_SCOOL="imputed_pad${PAD}_std${STD}_rp${RP}_r${R}_${SCOOL_BASENAME}.scool"

    python convertImputedToScool.py \
        --imputed-dir "$FILTERED_DIR" \
        --output-scool "$OUTPUT_SCOOL" \
        --chrom-sizes "$CHROM_SIZES_FILE" \
        --resolution "$RESOLUTION" \
        --chromosomes $CHROM_LIST

    if [ $? -eq 0 ]; then
        OUTPUT_SCOOLS+=("$OUTPUT_SCOOL")
        echo "  ✅ Created: $OUTPUT_SCOOL"
    else
        echo "  ❌ Failed to create scool for r = $R"
    fi
done

echo ""
echo "=========================================="
echo "Step 5: Verifying node degrees..."
echo "=========================================="

# Run node degree calculation for each generated scool
for SCOOL in "${OUTPUT_SCOOLS[@]}"; do
    echo ""
    echo "Testing: $(basename $SCOOL)"
    echo "----------------------------------------"
    python calc_node_degree.py --scool-file "$SCOOL"
done

echo ""
echo "=========================================="
echo "Imputation Complete!"
echo "=========================================="
echo "Generated ${#OUTPUT_SCOOLS[@]} scool file(s):"
for SCOOL in "${OUTPUT_SCOOLS[@]}"; do
    echo "  - $(basename $SCOOL)"
done
echo ""
echo "To use in the main pipeline:"
for SCOOL in "${OUTPUT_SCOOLS[@]}"; do
    echo "  ./example_run.sh imputation/$(basename $SCOOL)"
done
echo "=========================================="
