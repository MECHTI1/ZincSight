# Rotamer and execution-order diagnosis — 2026-10-30

Histidine rotamers were enabled in the notebook-default comparisons. The
notebook sets `include_histidine_rotamers = True`; diagnostic runs explicitly
set the equivalent engine option. The private preview entry point also enables
sampling by default. General MCP and repository CLI calls retain their opt-in
default: use `histidine_rotamers=true` or `--histidine-rotamers` when reproducing
the notebook. All counts below use the notebook's probability cutoff of 0.0001.

## Controlled native PostgreSQL runs

[Diagnostic run 41053002048](https://github.com/MECHTI1/ZincSight/actions/runs/41053002048)
uses six freshly downloaded public structures, two workers, Python 3.12 and the
pinned integration dependencies. It compares original scientific source
`967bdf431e3c86f0b390795022638eae32e64482` with three controlled variants:
only the local exclusion-set fix, only explicit ordering, and both changes.
Read-only tracing records candidate tables and actual worker inputs. Temporary
source copies leave the repository's scientific implementation untouched.

| Public input | Original, rotamers on | State fix only, on | Ordering only, on | Both, on | Original, off | Both, off |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1KLS | 1 | 1 | 1 | 1 | 1 | 1 |
| 2A0S-assembly1 | 33 | 37 | 33 | 36 | 14 | 12 |
| 8QEP | 27 | 27 | 27 | 27 | 18 | 18 |
| AF-P0A6G5-F1-model_v6 | 2 | 2 | 2 | 2 | 2 | 2 |
| AF-A0A068N621-F1-model_v6 | 5 | 5 | 5 | 5 | 5 | 5 |
| MGYP002718891411 | 2 | 2 | 2 | 2 | 2 | 1 |
| **Total** | **70** | **74** | **70** | **73** | **42** | **39** |

These are observed counts for this execution, not invariants of the unfixed
code. Unspecified row order and retained worker state can alter its output.
Equal counts also do not mean equal residue sets, coordinates or probabilities.

## Where results first diverge

All four rotamers-enabled native variants retain 1,271 initial motif candidates
and 494 after RMSD deduplication. Comparing original versus state-fix-only
records by structure, motif and ordered residues finds identical candidates and
RMSDs. The first substantive change is in rotamer geometry: 150 of the 494
records change a geometry/scoring field by more than 1e-6. Those changes feed
scoring and later proximity grouping.

The original `return_nonredundant_dicts_his_comb.py` stores excluded conformer
indices in a module-level set. An index can therefore remain excluded during a
later call, although that call supplies a different candidate list. Moving the
set inside `main` scopes exclusions to the list being compared. The geometry
formula and trained probability model are unchanged.

## Isolated history replay

`tests/replay_rotamer_history.py` exchanges only that helper while keeping the
geometry implementation and atom coordinates identical. For internal template
15 (DHH), it tests all 19 query candidates against 20 possible preceding inputs
captured from a fresh public-input trace, including the template input.

| Helper | Trials | Results changed by preceding input | Largest zinc-position change |
| --- | ---: | ---: | ---: |
| Original global exclusions | 380 | 187 | 1.36415475 Å |
| Exclusions local to the call | 380 | 0 | 0 Å |

There were no calculation errors in either mode. Change detection uses a 1e-8
threshold for zinc position or histidine-angle differences. This experiment
demonstrates history dependence; it does not measure biological accuracy.

A separate replay of the native run's 19 DHH candidates, restoring each recorded
initial exclusion set, reproduces its original geometry within 1.8e-9 Å in zinc
position and 2e-8 degrees in angles. The same check passes for the corrected
native run. Thus the captured state explains the freshly observed worker
behavior, rather than only a synthetic example.

Run the history check after collecting a trace with worker inputs:

```sh
python tests/replay_rotamer_history.py \
  --trace-dir /path/to/traces/current \
  --original-root /path/to/original-source \
  --output /path/to/history-replay.json
```

## Ordering and backend consistency

With the state fix applied, native original-order execution retains one extra
overlapping 2A0S prediction. Explicit ordering changes greedy proximity grouping
to retain 73 final sites. Its 73 shared sites have the same probabilities and
template IDs as the state-fix-only native run, with small numeric differences in
geometry. Stable grouping is required in addition to isolating rotamer state.

The combined native and PGlite traces retain the same 73 residue sets and
template IDs. Maximum probability difference is 2e-11 and maximum zinc-position
difference is 7e-9 Å, consistent with floating-point differences between hosts.

## Release limitation

The new checks explain why sampling being enabled does not guarantee stable
output in the original implementation. They establish a reproducibility defect
and verify the implemented fix on these inputs. They do not reconstruct an
unrecorded historical worker schedule, establish equivalence to every previous
result, or show that the revised search is biologically more accurate. Scientific
review and a benchmark against known sites remain necessary before treating the
preview as a replacement for the published scientific release.
