# ZincSight Preview

Explore potential zinc-binding sites in protein structures from a conversation.
Ask for a PDB ID, biological assembly, UniProt accession, AlphaFold model ID,
MGYP ID, or provide a PDB/mmCIF file. Results include coordinating residues,
predicted zinc coordinates, probability scores and template identifiers.

**Experimental integration:** outputs can differ from historical ZincSight
notebook results. Software checks establish reproducibility on tested inputs,
not the published method's benchmark accuracy for this revision.

## Try it

After installing the skill or plugin in a supported execution environment, ask:

> Use ZincSight Preview on PDB 1KLS. Report the coordinating residues, zinc
> position and probability, and save the full prediction JSON.

Histidine rotamer sampling is on by default, with the notebook's 0.0001 output
cutoff. A request processes one structure, up to 5 MB and 4,000 protein residues,
using the first structural model. This preview requires an explicit request.

The package contains one skill and its launcher. It starts embedded PostgreSQL
inside the task. It requires no ZincSight account, API key, external database
account or permanently running server. The host's plan, execution limits and
network permissions still apply; availability in every ChatGPT session is not
guaranteed.

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
public ChatGPT directory. Installation in one account does not make ZincSight
available in other accounts or cause automatic promotion.

Source and support: [MECHTI1/ZincSight](https://github.com/MECHTI1/ZincSight).
Report reproducible issues at [GitHub Issues](https://github.com/MECHTI1/ZincSight/issues),
including the software revision and effective settings, without private files.
No new license or service terms are imposed by this packaging change.
