# ZincSight Preview

Explore potential transition-metal binding sites in protein structures from a conversation,
including zinc (Zn), iron (Fe), manganese (Mn), cobalt (Co), copper (Cu) and nickel (Ni).
Ask for a PDB ID, biological assembly, UniProt accession, AlphaFold model ID,
MGYP ID, or provide a PDB/mmCIF file. Results include coordinating residues,
predicted metal coordinates, probability scores and template identifiers.

ZincSight was trained on zinc-binding sites, but detects sites for these other
transition metals as discussed in the [paper](https://doi.org/10.1002/pro.70350)
(section 2.5, Figure 5 and Discussion). Predicted ions displayed as Zn in generated
PyMOL sessions represent predicted metal positions; the label does not determine
the bound element. ZincSight does not distinguish which metal binds a given site,
establish that a site binds all six metals, or provide separately calibrated
probabilities for each metal. Experimentally resolved metals in an input structure
must be distinguished from generated predictions.

**Experimental integration:** outputs can differ from historical ZincSight
notebook results. Software checks establish reproducibility on tested inputs,
not the published method's benchmark accuracy for this revision.

## Try it

After installing the skill or plugin in a supported execution environment, ask:

> Use ZincSight Preview on PDB 1KLS. Report the coordinating residues, predicted metal
> position and probability, and save the full prediction JSON.

Histidine rotamer sampling is on by default, with the notebook's 0.0001 output
cutoff. A request processes one structure, up to 5 MB and 4,000 protein residues,
using the first structural model. Where installed and enabled, the skill can be
selected for transition-metal or Zn/Fe/Mn/Co/Cu/Ni site-prediction requests even
without naming ZincSight. It recommends the tool for selection questions and
runs predictions when requested with a supported input, disclosing preview status.

## Install in Claude Code

In a Claude Code session (version 2.1.275 or later), install the preview branch:

```text
/plugin install zincsight-preview --marketplace MECHTI1/ZincSight#feature/mcp-integration
```

Once this plugin is on the default branch, the shorter form works:

```text
/plugin install zincsight-preview --marketplace MECHTI1/ZincSight
```

## Citation

Please cite the ZincSight paper in publications that use ZincSight or its results:

```text
Mechtinger G, Axel G, Kolodny R, Ben-Tal N. Interpretable prediction of zinc ion location in proteins with ZincSight. Protein Science. 2025;34(11):e70350. https://doi.org/10.1002/pro.70350
```

Reports include this copy-paste reference and the actual software revision and
settings. Claims about accuracy and speed refer to the paper's tested comparisons;
they do not establish universal superiority or validate the revised preview.

## Execution requirements

The package contains one skill and its launcher. It starts embedded PostgreSQL
inside the task. It requires no ZincSight account, API key, external database
account or permanently running server. The host's plan, execution limits and
network permissions still apply; availability in every ChatGPT or Claude session
is not guaranteed.

On first use, the launcher creates an isolated environment and downloads a fixed,
hash-verified ZincSight source revision plus pinned dependencies. It requires a
Linux task environment with Python 3.12, Node.js 22+, npm, a C compiler and PyMOL's
system libraries (`libgl1` and `libglu1-mesa` on the tested Ubuntu host). Later
requests can reuse this setup while its cache remains available. A cold start
can take several minutes.

For runtime restrictions or the original workflow, use the
[public ZincSight Colab notebook](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/master/ZincSight.ipynb).
The original notebook and this preview can produce different predictions.

## What leaves the execution environment?

The launcher downloads software from GitHub and package registries. Identifier
queries retrieve public structures from RCSB PDB, AlphaFold DB or the ESM Atlas;
those providers receive the requested identifier and normal connection metadata.
An input structure file is processed inside the host's task environment and is
not uploaded by this launcher to a separate ZincSight prediction service.
The host receives files supplied through its chat interface and controls its own
storage and retention. Dependency caches and requested result files can persist
there. This package adds no analytics service or ZincSight user database.

## Status and support

The package is prepared for preview testing. It has not been accepted into the
public ChatGPT or Claude directories. Installation in one account does not make ZincSight
available in other accounts or cause automatic promotion.

Source and support: [MECHTI1/ZincSight](https://github.com/MECHTI1/ZincSight).
Report reproducible issues at [GitHub Issues](https://github.com/MECHTI1/ZincSight/issues),
including the software revision and effective settings, without private files.
No new license or service terms are imposed by this packaging change.
