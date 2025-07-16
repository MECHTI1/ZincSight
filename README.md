ZincSight
=========
**Interpretable prediction of zinc ion locations in proteins**

ZincSight is a tool for predicting zinc ion binding sites in protein structures. It is designed for **high-throughput analysis**, while keeping predictions **accurate** and **interpretable**.  

## Features

- Fast processing: Approximately 0.5–1 second per protein query on Google Colab—achieved by utilizing two threads on an Intel Xeon® CPU @ 2.20GHz, or alternatively, using two local 11th Gen Intel® Core™ i7-1165G7 CPU cores.
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

## 1. Quick Start with Google Colab
The fastest way to get started with ZincSight is through our interactive Google Colab notebook:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/master/ZincSight.ipynb)

## 2. Local Installation

### Downloading the repository
```bash
git clone https://github.com/MECHTI1/ZincSight.git
cd ZincSight
```

### Conda environment setup
```bash
conda env create -f environment.yml
conda activate zincsight
```

## 3. Running ZincSight Locally

After installing the environment, you can run ZincSight using the `zincsight.py` script with the following arguments:

| Argument                                          | Description                                  |
|---------------------------------------------------|----------------------------------------------|
| `-r`, `--rotamers`                                | Include histidine rotamers (default: no)     |
| `-i <IDENTIFIERs>`, `--identifiers <IDENTIFIERS>` | Comma-separated identifiers                  |
| `-f <FILE>`, `--file <FILE>`                      | Path of a text file with identifier per line |

At least one of `IDENTIFIERS` or `FILE` is required.

Script execution with the included example input file:
```shell
python ./zincsight.py -f example_input.txt
```
