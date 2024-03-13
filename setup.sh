#!/bin/bash

# Function to check if a command exists
command_exists () {
    type "$1" &> /dev/null ;
}

# Function to install MUSCLE
install_muscle() {
    echo "Installing MUSCLE..."
    wget https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
    tar -zxvf muscle3.8.31_i86linux64.tar.gz
    mv muscle3.8.31_i86linux64 muscle
    chmod +x muscle
    # Assuming ~/bin exists and is in your PATH
    mv muscle ~/bin
    echo "MUSCLE installed."
}

# Function to install BLAST
install_blast() {
    echo "Installing BLAST..."
    wget ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-aarch64-linux.tar.gz
    tar -zxvf ncbi-blast-2.15.0+-aarch64-linux.tar.gz
    mv ncbi-blast-2.15.0+ /opt/
    ln -s /opt/ncbi-blast-2.15.0+/bin/* ~/bin/
    echo "BLAST installed."
}

# Check if MUSCLE is installed
if command_exists muscle; then
    echo "MUSCLE is already installed."
else
    install_muscle
fi

# Check if BLAST is installed
if command_exists blastp; then
    echo "BLAST is already installed."
else
    install_blast
fi

# Define the name of the Conda environment
ENV_NAME="intron_conserv_env"

echo "Setting up the $ENV_NAME environment..."

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    echo "Conda could not be found. Please install Miniconda or Anaconda before continuing."
    exit 1
fi

# Create the Conda environment from the environment.yml file
conda env create -f environment.yml --name $ENV_NAME

echo "Please activate the $ENV_NAME environment using: conda activate $ENV_NAME"

echo "Setup completed."

