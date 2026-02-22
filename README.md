ZincSight
=========
**Interpretable prediction of zinc ion locations in proteins**

ZincSight is a tool for predicting zinc ion binding sites in protein structures. It is designed for **high-throughput analysis**, while keeping predictions **accurate** and **interpretable**.  

## Features

### Flexible input options
1. Protein IDs (e.g., PDB, AlphaFold, or ESM Metagenomic Atlas) 
2. Uploaded structure files in PDB or MMCIF format

### Comprehensive output:
  - Detailed CSV files including (among others):
    * predicted zinc-ion x,y,z coordinates
    * predicted binding residues
    * probability for being TP site
    * references for experimentally solved zinc-binding templates used in prediction
  - PyMOL sessions for easy visualization

## Quick Start with Google Colab
The fastest way to get started with ZincSight is through our interactive Google Colab notebooks:

- **ZincSight Colab:** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/master/ZincSight.ipynb)
- 
- **ZincSight MaxDrive Colab:** Built for massively scalable screening across millions of structures, with native Google Drive integration for large input/output workflows  
  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/master/ZincSight_MaxDrive.ipynb)

  
## Local Installation

### Downloading the repository
```bash
git clone https://github.com/MECHTI1/ZincSight.git
cd ZincSight
```

### Dependencies installation
```bash
pip install -r requirements.txt
```

### ***Optional***: adding a new zinc-binding motif template to the training set

In order to add a new zinc-binding motif from a PDB structure that was not available at the time of training, use the `add_template.py` script with the following arguments:

| Argument | Description                             |
|----------|-----------------------------------------|
| PDB_ID   | PDB identifier                          |
| RESIDUES | Comma-separated pairs of CHAIN_RESIDUEs |

Example execution:
```shell
python ./add_templates 4f3w B_89,B_56,B_92
```

### ZincSight execution

After installing the environment, you can run ZincSight using the `zincsight.py` script with the following arguments:

| Argument                                          | Description                                  |
|---------------------------------------------------|----------------------------------------------|
| `-r`, `--rotamers`                                | Include histidine rotamers (default: no)     |
| `-i <IDENTIFIERS>`, `--identifiers <IDENTIFIERS>` | Comma-separated identifiers                  |
| `-f <FILE>`, `--file <FILE>`                      | Path of a text file with identifier per line |
| `-c <NUMBER>`, `--cores <NUMBER>`                 | Number of cores (default: 2)                 |


At least one of `IDENTIFIERS` or `FILE` is required.

Script execution with the included example input file:
```shell
python ./zincsight.py -f example_input.txt
```
