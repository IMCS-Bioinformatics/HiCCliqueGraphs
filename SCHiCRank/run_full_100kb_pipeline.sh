#!/bin/bash

# ============================================================================
# run_full_100kb_pipeline.sh
# Full pipeline: Resolution change -> Imputation -> Motif ordering
# ============================================================================

set -e

# Configuration
SOURCE_SCOOL="sourceData/nagano_10kb_cell_types.scool"
CELL_TYPE_FILE="sourceData/nagano_assoziated_cell_types.txt"
COARSEN_FACTOR=10
IMPUTATION_R=3.0

# Derived paths
COARSENED_SCOOL="resolutionChange/coarsened_k${COARSEN_FACTOR}_nagano_10kb_cell_types.scool"
IMPUTED_SCOOL_NAME="imputed_pad1_std1_rp0.5_r${IMPUTATION_R}_coarsened_k${COARSEN_FACTOR}_nagano_10kb_cell_types.scool"
IMPUTED_SCOOL="imputation/${IMPUTED_SCOOL_NAME}"

# Motif ordering config (imputed file is already at 100kb, so k=1)
POSTFIX="imputed100k_r${IMPUTATION_R}"
K_MULTIPLIER=1
FN_RESOLUTION=100000
RESOLUTION=100000

WORK_DIR="motifOrdering/${POSTFIX}"
SCRIPT_DIR="motifOrdering"

# Activate conda environment
echo "Activating schicluster environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
CONDA_ENV_NAME=$(grep '^CONDA_ENV_NAME=' setup_schicluster_env.sh | cut -d'"' -f2)
conda activate "$CONDA_ENV_NAME"

echo ""
echo "================================================================"
echo "Full 100kb Pipeline: Coarsen -> Impute -> Motif Ordering"
echo "================================================================"
echo "Source:         $SOURCE_SCOOL"
echo "Coarsen factor: $COARSEN_FACTOR (10kb -> 100kb)"
echo "Imputation r:   $IMPUTATION_R"
echo "Postfix:        $POSTFIX"
echo "Work dir:       $WORK_DIR"
echo "================================================================"
echo ""

# ==========================================
# Step 1: Resolution change (10kb -> 100kb)
# ==========================================
echo "========================================"
echo "Step 1: Resolution change to 100kb"
echo "========================================"

if [ -f "$COARSENED_SCOOL" ]; then
    echo "Coarsened file already exists: $COARSENED_SCOOL"
    echo "Skipping resolution change."
else
    cd resolutionChange
    python coarsen_scool.py \
        --input-scool "../$SOURCE_SCOOL" \
        --output-scool "coarsened_k${COARSEN_FACTOR}_nagano_10kb_cell_types.scool" \
        --factor "$COARSEN_FACTOR"
    cd ..
fi

echo ""

# ==========================================
# Step 2: Imputation
# ==========================================
echo "========================================"
echo "Step 2: Imputation with r=$IMPUTATION_R"
echo "========================================"

if [ -f "$IMPUTED_SCOOL" ]; then
    echo "Imputed file already exists: $IMPUTED_SCOOL"
    echo "Skipping imputation."
else
    cd imputation
    ./run_full_imputation.sh "../$COARSENED_SCOOL" ALL "$IMPUTATION_R"
    cd ..
fi

# Check imputed file exists
if [ ! -f "$IMPUTED_SCOOL" ]; then
    echo "Error: Imputed file not found: $IMPUTED_SCOOL"
    echo "Check the imputation output for errors."
    exit 1
fi

echo ""

# ==========================================
# Step 3: Motif ordering
# ==========================================
echo "========================================"
echo "Step 3: Motif ordering"
echo "========================================"
echo ""

mkdir -p "$WORK_DIR"

# Step 3.0: Generate cell phase metadata
python "$SCRIPT_DIR/generateCellPhaseInfo.py" \
    --scool-file "$IMPUTED_SCOOL" \
    --cell-type-file "$CELL_TYPE_FILE" \
    --output cellAndPhaseInfo.pkl \
    --work-dir "$WORK_DIR"

# Step 3.1: Process scool into per-chromosome interaction files
python "$SCRIPT_DIR/runProcessCoolDataset.py" \
    -k "$K_MULTIPLIER" \
    --chromosomes ALL \
    --fn "$IMPUTED_SCOOL" \
    --fn-resolution "$FN_RESOLUTION" \
    --postfix "$POSTFIX" \
    --work-dir "$WORK_DIR"

# Step 3.2: Generate clique data (K3-K8)
python "$SCRIPT_DIR/runCreateCliques.py" \
    --chromosomes ALL \
    --postfix "$POSTFIX" \
    --resolution "$RESOLUTION" \
    --work-dir "$WORK_DIR"

# Step 3.3: Compute pairwise similarities
python "$SCRIPT_DIR/runCreatePairwiseSimilarities.py" \
    --chromosomes ALL \
    --postfix "$POSTFIX" \
    --resolution "$RESOLUTION" \
    --work-dir "$WORK_DIR"

K_NEIGHBORS=5
MIN_ACTIVE_CELLS=10
MIN_MOTIF_COUNT=10

# PageRank method — cliques
echo ""
echo "========================================"
echo "Computing PageRank orderings (cliques)..."
echo "========================================"

for K in K3 K4 K5 K6 K7 K8; do
    for L in alllengths long; do
        MOTIF_CONFIG="${K}-${L}"
        echo ""
        echo "Processing motif: $MOTIF_CONFIG"
        python "$SCRIPT_DIR/runPageRankFilter.py" \
            --input-dir "${WORK_DIR}/pairwiseSimilarities/${MOTIF_CONFIG}/" \
            --label "pairwiseSimilarities_${MOTIF_CONFIG}" \
            --k "$K_NEIGHBORS" \
            --min-active-cells "$MIN_ACTIVE_CELLS" \
            --no-plots \
            --work-dir "$WORK_DIR"
        python "$SCRIPT_DIR/runCalculateTau.py" \
            --input-csv "final_active_cells_pairwiseSimilarities_${MOTIF_CONFIG}.csv" \
            --motif-config "$MOTIF_CONFIG" \
            --k "$K_NEIGHBORS" \
            --min-active-cells "$MIN_ACTIVE_CELLS" \
            --log-file tau_scores.txt \
            --work-dir "$WORK_DIR"
    done
done

# PageRank method — interactions
echo ""
echo "========================================"
echo "Computing PageRank orderings (interactions)..."
echo "========================================"

python "$SCRIPT_DIR/createInteractionPairwiseSimilarities.py" \
    --work-dir "$WORK_DIR" \
    --config "$POSTFIX"
for L in alllengths long; do
    MOTIF_CONFIG="interaction-${L}"
    echo ""
    echo "Processing motif: $MOTIF_CONFIG"
    python "$SCRIPT_DIR/runPageRankFilter.py" \
        --input-dir "${WORK_DIR}/pairwiseSimilarities/${MOTIF_CONFIG}/" \
        --label "pairwiseSimilarities_${MOTIF_CONFIG}" \
        --k "$K_NEIGHBORS" \
        --min-active-cells "$MIN_ACTIVE_CELLS" \
        --no-plots \
        --work-dir "$WORK_DIR"
    python "$SCRIPT_DIR/runCalculateTau.py" \
        --input-csv "final_active_cells_pairwiseSimilarities_${MOTIF_CONFIG}.csv" \
        --motif-config "$MOTIF_CONFIG" \
        --k "$K_NEIGHBORS" \
        --min-active-cells "$MIN_ACTIVE_CELLS" \
        --log-file tau_scores.txt \
        --work-dir "$WORK_DIR"
done

echo ""
echo "========================================"
echo "Computing motif-count-based orderings..."
echo "========================================"

for K in K3 K4 K5 K6 K7 K8; do
    for L in alllengths long; do
        MOTIF_CONFIG="${K}-${L}"

        echo "Processing motif-count ordering: $MOTIF_CONFIG"

        python "$SCRIPT_DIR/calculateTauFromMotifCounts.py" \
            --motif-type "$MOTIF_CONFIG" \
            --postfix "$POSTFIX" \
            --resolution "$RESOLUTION" \
            --min-count "$MIN_MOTIF_COUNT" \
            --log-file tau_scores.txt \
            --work-dir "$WORK_DIR"
    done
done

echo ""
echo "========================================"
echo "Generating tau scores visualization..."
echo "========================================"

python "$SCRIPT_DIR/plotTauScores.py" \
    --input tau_scores.txt \
    --output tau_comparison.png \
    --work-dir "$WORK_DIR"

echo ""
echo "================================================================"
echo "Pipeline complete!"
echo "================================================================"
echo "Results in: $WORK_DIR/"
echo "  - tau_scores.txt        (all tau scores)"
echo "  - tau_comparison.png    (visualization)"
echo "================================================================"
