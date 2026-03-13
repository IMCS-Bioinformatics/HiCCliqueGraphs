#!/bin/bash

# ScHiCluster environment setup script
# Creates the conda environment if it doesn't exist, then activates it

CONDA_ENV_NAME="schicrank_env"

# Initialize conda
source "$(conda info --base)/etc/profile.d/conda.sh"

# Check if the environment exists
if ! conda env list | grep -q "^${CONDA_ENV_NAME} "; then
    echo "Creating conda environment: $CONDA_ENV_NAME"
    conda create -n "$CONDA_ENV_NAME" python=3.6.8 r-base -c conda-forge -y
else
    echo "Environment $CONDA_ENV_NAME already exists."
fi

echo "Activating environment: $CONDA_ENV_NAME"
conda activate "$CONDA_ENV_NAME"

# Check if scHiCluster is installed
if ! pip show scHiCluster &> /dev/null; then
    echo "Upgrading pip and setuptools..."
    pip install --upgrade pip setuptools wheel

    echo "Installing scHiCluster and dependencies..."
    # Force rpy2 to use ABI mode to avoid needing R installed
    # export RPY2_CFFI_MODE=ABI
    pip install git+https://github.com/zhoujt1994/scHiCluster.git

    echo "Installing opencv-python (compatible version for Python 3.7)..."
    pip install opencv-python==4.5.5.64

    echo "Installing additional dependencies..."
    pip install h5py pandas numpy scipy networkx matplotlib

    echo "Environment setup complete!"
else
    echo "scHiCluster already installed."
fi

# Verify installation
echo ""
echo "Testing scHiCluster installation..."
hicluster --help
