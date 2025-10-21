#!/bin/bash

# Setup script for TL1A PET Imaging Agent project

echo "Setting up Python environment..."

# Create conda environment
conda create -n tl1a_env python=3.10 -y

# Activate environment
conda activate tl1a_env

# Install dependencies
pip install -r requirements.txt

# Download BioTransformer JAR (placeholder - user needs to download manually)
echo "Please download BioTransformer3.0.jar from https://github.com/BioTransformer/BioTransformer3.0-cli/releases and place it in scripts/"

echo "Setup complete. Activate with: conda activate tl1a_env"