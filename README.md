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
- **ZincSight BatchDrive Colab:** Built for massively scalable screening across millions of structures, with native Google Drive integration for large input/output workflows  
  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/master/ZincSight_MaxDrive.ipynb)

  
## Local Installation

### Downloading the repository
```bash
git clone https://github.com/MECHTI1/ZincSight.git
cd ZincSight
```

### Dependencies installation
Tested with Python 3.12 and 3.13 on Linux:
```bash
pip install -r requirements.txt
```
Downloading AlphaFold or TED models by ID also needs `aria2c` (`apt-get install aria2` or `brew install aria2`).

### ***Optional***: adding a new zinc-binding motif template to the training set

In order to add a new zinc-binding motif from a PDB structure that was not available at the time of training, use the `add_template.py` script with the following arguments:

| Argument | Description                             |
|----------|-----------------------------------------|
| PDB_ID   | PDB identifier                          |
| RESIDUES | Comma-separated pairs of CHAIN_RESIDUEs |

Example execution:
```shell
python add_template.py 4f3w B_89,B_56,B_92
```

### ZincSight execution

After installing the environment, you can run ZincSight using the `zincsight.py` script with the following arguments:

| Argument                                          | Description                                  |
|---------------------------------------------------|----------------------------------------------|
| `-r`, `--rotamers`                                | Include histidine rotamers (default: no; the Colab notebooks use them) |
| `-i <IDENTIFIERS>`, `--identifiers <IDENTIFIERS>` | Comma-separated identifiers                  |
| `-f <FILE>`, `--file <FILE>`                      | Path of a text file with identifier per line |
| `-s <PATH> ...`, `--structures <PATH> ...`         | Local PDB/mmCIF files, or folders of them    |
| `-c <NUMBER>`, `--cores <NUMBER>`                 | Number of cores (default: 2)                 |
| `-p`, `--pse_output`                              | Also write PyMOL `.pse` sessions             |

At least one of `IDENTIFIERS`, `FILE` or `STRUCTURES` is required. ZincSight predicts on every file in its `Query AlphaFold structures` folder, so empty that folder between unrelated runs. Settings such as the probability cutoff can be overridden through environment variables; see [.env.example](.env.example).

Script execution with the included example input file:
```shell
python ./zincsight.py -f example_input.txt
```

## AI assistant integration (development preview)

**ZincSight Preview** now has a one-skill package for on-demand predictions in a
supported cloud task, with embedded PostgreSQL and no separate API key or hosted
database. See [installation and sharing](integrations/mcp/DISTRIBUTION.md) for the
repository marketplace, tested package and example prompts. The same package
installs in Codex and Claude. This experimental package is not yet listed in the
public ChatGPT or Claude directories.

An optional MCP wrapper for ChatGPT, Claude and other compatible clients is
available in [integrations/mcp](integrations/mcp/README.md). It includes isolated
request processing and structured prediction results. It is not yet a hosted
service or a published plugin; see the integration guide for validation status
and public-launch requirements.
The optional embedded PostgreSQL mode also supports direct predictions inside
a cloud task, without permanent hosting; see the guide for setup and limitations.

## Citation and license

If you use ZincSight or its results, please cite:

> Mechtinger G, Axel G, Kolodny R, Ben-Tal N. Interpretable prediction of zinc ion location in proteins with ZincSight. *Protein Science.* 2025;34(11):e70350. https://doi.org/10.1002/pro.70350

GitHub's "Cite this repository" button gives the same reference from [CITATION.cff](CITATION.cff).

ZincSight is free to use under the [Creative Commons Attribution-NoDerivatives 4.0](LICENSE) license (CC BY-ND 4.0). You may use and share it unmodified for any purpose. You may not distribute modified versions.
