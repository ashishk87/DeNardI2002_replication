# DeNardi_2004_replication

This repository replicates the De Nard and Tristani (2002) study using Julia. It includes all necessary scripts and a Jupyter Notebook for analysis.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running the Notebook](#running-the-notebook)

## Prerequisites

Ensure the following are installed on your local machine:

- [Julia 1.10.2](https://julialang.org/downloads/)
- [Git](https://git-scm.com/downloads)
- [JupyterLab](https://jupyter.org/install)

## Installation

Follow these steps to set up the environment and run the project locally.

### 1. Clone the Repository

```bash
git clone https://github.com/ashishk87/DeNardI2002_replication
cd your-repository
```

### 2. Run the Setup Script
A ```setup.sh``` script is provided to automate the environment setup. Run the following command:
```
./setup.sh
```

If you encounter a "Permission denied" error, make the script executable:
```
chmod +x setup.sh
./setup.sh
```

### Running the Notebook
In JupyterLab (or VSCode), navigate to the src folder and open DeNardi200_Replication.ipynb. Run the cells to execute the replication study.






