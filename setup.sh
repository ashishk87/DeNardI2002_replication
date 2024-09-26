#!/bin/bash

# Ensure the script exits if any command fails
set -e

# Check if Julia is installed
if ! command -v julia &> /dev/null
then
    echo "Julia could not be found. Please install Julia from https://julialang.org/downloads/"
    exit
fi

# Activate the Julia environment and instantiate
echo "Activating Julia environment..."
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Check if IJulia is installed
echo "Checking for IJulia installation..."
julia --project=. -e '
    if !("IJulia" in names(Pkg.installed()))
        println("IJulia not found. Installing IJulia...")
        Pkg.add("IJulia")
    else
        println("IJulia is already installed.")
    end
'

# Build IJulia
echo "Building IJulia..."
julia --project=. -e 'using Pkg; Pkg.build("IJulia")'

echo "Setup complete! You can now run JupyterLab with 'jupyter lab'."
