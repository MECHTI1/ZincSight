# Notebook execution and consistency checks — 2026-10-30

**Status: development preview. Runtime consistency is tested; scientific release
review is still required.**

The subsequent [rotamer diagnosis](ROTAMER_DIAGNOSIS.md) confirms that sampling
was enabled, isolates retained exclusion state using 380 controlled history
trials, and records native PostgreSQL comparisons with sampling on and off.

## Public notebook configuration

`ZincSight.ipynb` clones this repository. The tracked `.env` enables
`FILTER_PROB_AFTER_COMPRESSION=true` and sets
`EXCLUDE_PREDICTIONS_WITH_PROB_THRESHOLD=0.0001`. Histidine sampling is enabled
in the notebook's example input cell.

The initial reproduction scripts disabled `.env` and probability filtering.
Their totals of 137 predictions for the revised integration and 138 for the
original source were therefore **unfiltered engine totals**, not the notebook's
default output. Applying the public notebook's cutoff to those saved public-input
runs produces:

| Public input | Original native code | Current corrected integration |
| --- | ---: | ---: |
| 1KLS | 1 | 1 |
| 2A0S-assembly1 | 33 | 36 |
| 8QEP | 27 | 27 |
| AF-P0A6G5-F1-model_v6 | 2 | 2 |
| AF-A0A068N621-F1-model_v6 | 5 | 5 |
| MGYP002718891411 | 2 | 2 |
| **Total** | **70** | **73** |

These figures derive from newly computed public structures. Equal counts alone
do not imply equal predictions or establish biological accuracy.

The wrapper, MCP schema and direct CLI now default to the notebook's 0.0001
cutoff. Set `minimum_probability=0` or CLI `--minimum-probability 0` to inspect
all computed sites. Filtering occurs after prediction and before truncation;
each result records the effective cutoff. This is an output filter, not a
biological confidence claim. General CLI/MCP histidine sampling remains opt-in;
the private preview skill enables it by default.

## Implemented reproducibility changes

- Explicit SQL ordering for candidate/residue reads, template traversal, RMSD
  coordinate pairs and greedy proximity grouping stabilizes unspecified row order.
- The histidine-rotamer exclusion set is now local to each comparison. Previously
  excluded indices could leak into later calls, remove unrelated conformers or
  cause an IndexError. Three focused tests fail before this fix and pass afterward.
- First-model mmCIF preparation preserves original coordinate text, avoiding
  rounding to three decimal places during serialization.
- Empty prediction output is handled cleanly, probability-calculation failures
  raise errors, and temporary motif tables are explicitly released for the
  embedded database's single-session adapter.

Geometry/scoring formulas and model weights are unchanged, but candidate
selection and state/order corrections can change predictions. They require
scientific review before a public release; consistency is not a new accuracy
benchmark.

## Completed validation before the cutoff correction

[Run 41036947949](https://github.com/MECHTI1/ZincSight/actions/runs/41036947949)
passed 45 tests on Ubuntu 24.04, Python 3.12 and Node.js 22. The six-example
workflow compared all 137 unfiltered predictions between native PostgreSQL and
PGlite, and between batch execution and individual wrapper calls. Coordinating
residue sets, probabilities, coordinates and template IDs agreed.

The unchanged scientific source at
`967bdf431e3c86f0b390795022638eae32e64482` was also run on native PostgreSQL with:

1. Python 3.12 and the integration's pinned dependencies.
2. Python 3.13 and the public notebook's dependency specification and pgserver/QCP
   compatibility shims, in
   [run 41038575388](https://github.com/MECHTI1/ZincSight/actions/runs/41038575388).

Both original-code runs had 138 unfiltered rows and 70 after applying 0.0001.
Their runtime versions and fresh prediction outputs are retained in workflow
artifacts. The original-code workflow now explicitly applies and records the
public notebook's filter settings during execution.

## Regression workflow

`tests/notebook_regression.py` downloads the six public examples from fixed
RCSB, AlphaFold and ESM providers. It records source URLs and input hashes, runs
isolated original/current native and current embedded jobs, and retains full
engine CSVs for backend comparisons. Reference and default-wrapper comparisons
apply the notebook's 0.0001 cutoff consistently. Differences in current backend
or batch/wrapper results fail the consistency check.

An optional `--expected /local/reference.csv` supports private local comparison;
that input is never uploaded by this script. CI uses only fresh public inputs
and newly computed outputs. Repository source, environment settings, input
coordinates, execution order and effective thresholds must all be recorded when
investigating scientific differences.

The 5 MB / 4,000-residue wrapper limit includes the 8QEP example, and the
180-second per-prediction limit remains. Identifier resolution includes PDB
assemblies, UniProt-to-AlphaFold v6, explicit AlphaFold model versions and MGYP.
The tested cloud runtime uses embedded PostgreSQL because native pgserver cannot
create its operating-system user there.

## Worker-count regression

The current embedded batch is also run with one worker. All unfiltered residue
sets, template IDs, probabilities and coordinates must match the two-worker
batch under the same comparison tolerances. On 1 November 2026, all 137 matched
exactly (73 after the notebook cutoff). The diagnostic workflow additionally
runs the small resolved-metal coordinate check described in
[the rotamer diagnosis](ROTAMER_DIAGNOSIS.md#additional-checks--2026-11-01).
