# Notebook-default regression check — 2026-10-30

**Status: historical-output equivalence has not passed. Keep this integration a
development preview until the differences are explained and resolved.**

**Follow-up:** the rotamer-state bug described below has now been fixed. All six
public examples agree across native PostgreSQL, PGlite and individual wrapper
calls. The wrapper also supports the notebook identifiers and 8QEP's size, and
preserves mmCIF coordinate precision. These fixes do not establish agreement with
the historical archive. See the follow-up results at the end of this report.

The notebook's six default query IDs are a useful existing reference set. The
maintainer supplied `ZincSight_sample_id_r-sampled_2026-10-14_08-54-27.tar.gz`,
containing 72 prediction rows, six PyMOL sessions and a template-reference CSV.
The notebook defaults enable histidine rotamers. The archive SHA-256 is
`7803673b27ec0e8b8dc66c79a5607cff072f607638512f2e6662cf132894b396`.

## Method

This check used the code at `579ecbe87b919908c036ae909957e7eb503ef1e8`,
Python 3.12, the integration's pinned Python dependencies, Node.js 24 and the
locked PGlite dependencies. Original molecule objects were recovered from the
uploaded sessions, excluding predicted-zinc objects and session commands, and
their first states exported as mmCIF. These are reconstructed structure inputs,
not byte-identical original downloads; export precision and alternate-location
handling need checking when interpreting small differences.

Histidine sampling was enabled, probability filtering disabled and all sites
retained. Site matching used the unordered set of coordinating chain/residue
identifiers; no duplicate residue sets occurred in these comparisons. Probabilities,
coordinates, template IDs and template RMSDs were compared for matching sites.

Five inputs ran through the actual embedded wrapper. 8QEP has 3,177 protein
residues and exceeds its 2,000-residue cap, so it was tested through the embedded
engine separately, with the same 180-second timeout. That engine test completed
in approximately 26 seconds; it does not establish wrapper support for 8QEP.
This describes the initial check; the limit has since been increased and 8QEP
validated through the wrapper as detailed below.

## Results

| Structure | Archived sites | Current single-input sites | Matching residue sets |
| --- | ---: | ---: | ---: |
| 1KLS | 1 | 1 | 1 |
| 2A0S-assembly1 | 35 | 60 | 29 |
| 8QEP (engine only) | 27 | 55 | 26 |
| AF-A0A068N621-F1-model_v6 | 5 | 14 | 5 |
| AF-P0A6G5-F1-model_v6 | 2 | 3 | 2 |
| MGYP002718891411 | 2 | 4 | 2 |

- 1KLS reproduces the archived probability and zinc position exactly, with only
  floating-point-scale RMSD difference.
- For 2A0S-assembly1, sites with probability above 0.9 decrease from two to zero.
  A shared site's probability changes from 0.55065227 to 0.03764699. This is a
  substantive difference, not a small numeric-tolerance issue.
- AF-A0A068N621-F1-model_v6 retains all five archived residue sets but adds nine
  sites. One shared probability changes from 0.040473867 to 0.31087896.
- AF-P0A6G5-F1-model_v6 and MGYP002718891411 reproduce probabilities and positions
  of their archived sites, but return additional low-probability sites.
- 8QEP retains five sites above probability 0.9, but its complete output differs.

The six structures were also tested together through the current embedded
engine. Counts remained 1, 60, 55, 14, 3 and 4 respectively, but some probabilities
changed from the single-input runs. In particular, the maximum archived-versus-
current probability difference for AF-A0A068N621-F1-model_v6 decreased from
0.270405093 to 0.002319346. A separate direct-engine run using the same named,
unmodified mmCIF file alone reproduced the single-input result, confirming that
the batch difference was not caused by the wrapper's renaming or normalization.

As a diagnostic only, the six ordering-modified source files were restored from
original `master` commit `967bdf431e3c86f0b390795022638eae32e64482` in a temporary
copy, retaining explicit temporary-table cleanup required by PGlite. That batch
returned 1, 65, 55, 15, 3 and 4 sites respectively. Simply restoring the former
SQL ordering therefore did not reproduce the archive either. This was still a
PGlite run, not a reproduction of the original notebook's native database and
dependency environment.

## Confirmed diagnostic issue and next checks

`return_nonredundant_dicts_his_comb.py` keeps
`processed_indices_to_move_out_dict_1` at module scope. A minimal two-call probe
confirmed that marking index zero as redundant in one call causes a distinct
conformer at index zero in the next call to be removed too; clearing the state
before that second call preserves it. This code predates the integration.
Its contribution to the full prediction discrepancies has not yet been isolated,
and the helper was not changed during that initial check.

The initial follow-up plan was:

1. Reproduce the reference run with original native PostgreSQL code and pinned
   inputs/settings, checking reconstructed coordinates against original downloads.
2. Isolate order, batch and rotamer-state effects; establish the intended behavior
   with the maintainer before accepting changes to historical predictions.
3. Add regression coverage for the resolved behavior, including single versus
   batch execution and repeated rotamer calculations.
4. Validate support for large notebook examples and the additional identifier
   sources. At that point, the tool directly resolved only four-character PDB
   IDs; other identifiers required file contents.

The earlier 31 passing integration tests remain useful but compare the revised
backends on other controls. They do not supersede these historical-output
differences or establish scientific equivalence. Neither counts nor agreement
with this archive alone establish biological prediction accuracy.

## Follow-up: consistent execution and original-notebook reproduction

The code tested at `b9f1b72a2a40cb6ac9c90aef6bf9ae6293c1c3da`:

- Keeps the rotamer exclusion set local to each call. Three targeted tests fail
  before the fix and pass afterward, including a previously possible IndexError.
- Preserves original mmCIF coordinate text while selecting the first model.
  This avoids rounding high-precision inputs to three decimal places.
- Resolves all notebook identifier types and accepts 5 MB / 4,000 residues, with
  the 180-second prediction timeout unchanged. 8QEP now passes through the actual
  wrapper, and the UniProt shorthand is exercised end to end.

GitHub Actions run [41036947949](https://github.com/MECHTI1/ZincSight/actions/runs/41036947949)
passed 45 tests plus the six-example workflow. Native and embedded databases
return identical sites/probabilities/coordinates across all 137 current
predictions. Each individual wrapper call also matches its full-batch result.
These comparisons use fresh public downloads, not a published copy of the
maintainer's archive. The archived CSV remains a local reference.

Two original-code reproductions were also performed on the untouched scientific
source at `967bdf431e3c86f0b390795022638eae32e64482`:

1. Python 3.12 with the integration's dependency set and native pgserver.
2. Python 3.13 with the notebook's dependency specification and its pgserver/QCP
   compatibility shims. This resolved to pixeltable-pgserver 0.7.1, Biopython 1.88,
   scikit-learn 1.9.1 and SciPy 1.18.1. The complete installed package list and
   predictions are preserved in the workflow artifact for
   [run 41038575388](https://github.com/MECHTI1/ZincSight/actions/runs/41038575388).

Both original-code runs produced these counts:

| Public input | Original code (both environments) | Current corrected integration |
| --- | ---: | ---: |
| 1KLS | 1 | 1 |
| 2A0S-assembly1 | 59 | 60 |
| 8QEP | 56 | 55 |
| AF-P0A6G5-F1-model_v6 | 3 | 3 |
| AF-A0A068N621-F1-model_v6 | 15 | 14 |
| MGYP002718891411 | 4 | 4 |

Thus the October archive mismatch cannot be attributed solely to the integration.
The current notebook setup does not reconstruct the exact historical run, and
identical counts do not imply identical site details. The integration also has
real changes relative to freshly run original code, consistent with its state
and ordering fixes. Do not claim historical equivalence or unchanged published
benchmark accuracy.

Fresh-download protein atom identities and coordinates match the recovered
session structures for five examples. 8QEP's session export introduced rounding
of up to approximately 0.00085 angstrom. Replacing the QCP calculation with the
notebook's SciPy compatibility implementation in an isolated diagnostic run did
not explain the major discrepancies.

Remaining scientific work is to identify the exact code/setup that generated the
October archive and agree on the intended regression behavior before public
release. The reproducibility fixes and runtime consistency checks are complete
on this six-example set; they are not a new biological accuracy benchmark.
