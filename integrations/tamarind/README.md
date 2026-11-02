# ZincSight on Tamarind

This branch packages **ZincSight Preview** as a Tamarind Custom Tool. It accepts
one uploaded PDB/mmCIF structure per job and writes `sites.csv` and
`prediction.json`. It runs on two CPU cores with embedded PostgreSQL (PGlite)
over a private, job-local Unix socket;
no GPU or external database is required.

## Which branch to connect

| Tamarind setting | Value |
| --- | --- |
| Repository | `MECHTI1/ZincSight` |
| Branch | `feature/tamarind-integration` |
| Source directory | Repository root |
| Suggested tool ID | `zincsight-preview` |
| Display name | `ZincSight Preview` |
| Automatic publishing | Off during initial validation |

The branch is based on `feature/mcp-integration`. It adds Tamarind's root-level
`Dockerfile`, `config.json`, and `run.sh`, plus this adapter and tests. The existing
`master`, `sqlite`, and MCP branches do not require edits to use this branch.
The scientific engine in the image is pinned to commit
`7bf80c27fd7f0e8ce7b6b100e305a7bc18620a03`, with its source archive SHA-256 verified
during the build. Changing the branch's scientific source alone does not change
the deployed engine: update the pin and digest in `Dockerfile` and the matching
`ENGINE_REVISION` in `predict.py`, then rerun validation.

The Docker image applies a checked transport-only patch to the pinned PGlite
launcher: its database socket lives inside the prediction's temporary directory,
so the worker does not depend on loopback TCP being reachable in Tamarind's
isolated runtime. The database fixture, scientific code, and trained weights
remain pinned. Results record `database_transport: "unix_socket"` and integration
version `0.1.2`. Failed jobs include bounded worker/database diagnostics in their
job log before temporary files are removed.

## Connect, build, test, publish

1. Open [Tamarind Custom Tools](https://app.tamarind.bio/app/custom-tools) in your
   account and create a tool with the suggested ID above (or another available ID).
2. In its **Settings** tab, link this repository and branch. An organization
   administrator may need to install Tamarind's GitHub App and grant access to
   this repository. Keep automatic publishing off for the first build.
3. Build/deploy the source. The root `config.json` defines the inputs and CPU
   resources. Dependencies and the pinned engine are installed under `/opt`,
   because Tamarind mounts submitted source over `/app` at runtime.
4. In **Test**, upload a small public structure such as
   [1KLS.cif](https://files.rcsb.org/download/1KLS.cif), enable histidine rotamers,
   and use the default cutoff `0.0001`. Confirm the build's test job finishes and
   its CSV/JSON files are available. A Tamarind test job uses account compute.
5. After reviewing that test, publish the completed version to your organization.
   Use **Manage Visibility** to choose which organization members can run it.

These are the current [Tamarind Custom Tool conventions](https://docs.tamarind.bio/tamarind/custom-tools).
Creating a GitHub branch does not itself create or publish a Tamarind tool.
An account with Custom Tools access is required. Custom tools are private to the
organization by default; public catalog inclusion requires a separate arrangement
with Tamarind. They can be used through Tamarind's UI, API, and MCP once deployed
and available to the caller.

## Inputs and results

- `structure`: one `.pdb`, `.cif`, or `.mmcif` upload. Maximum 5 MB and 4,000
  protein residues; the first structural model and all its protein chains are used.
- `histidine_rotamers`: defaults to `true`.
- `minimum_probability`: defaults to the notebook output cutoff, `0.0001`.
  This is an output filter, not a biological confidence threshold.
- `max_sites`: integer from 1 to 100, default 100. The highest-probability sites
  are returned. JSON includes `site_count` (total passing the cutoff),
  `returned_site_count`, and `truncated`; the CSV contains only returned sites.

Tamarind jobs have no internet access. Download structures before uploading;
PDB/UniProt IDs and external URLs are not accepted by this adapter. Results record
the input SHA-256, source revision, settings, model policy, limitations, and paper
citation. Empty results produce valid JSON and a header-only CSV. Invalid inputs
or engine failures produce a failed job rather than a successful empty prediction.

This adapter exports tabular and structured results; it does not export the
original notebook's PyMOL sessions. Coordinates are in the input structure's frame.

## Validation and scientific status

The **Tamarind validation** GitHub Actions workflow validates `config.json` against
Tamarind's live schema, builds the Docker image, and runs it with
`--network none`. It mounts only the adapter and inputs at `/app`, exercising the
same source-overlay behavior as Tamarind. Checks cover PDB/mmCIF, rotamers on/off,
truncation, an empty result, invalid input, and exact agreement with the pinned
wrapper for the same inputs and settings. This does not replace a test job on
Tamarind itself.

```bash
# Run from this branch's repository root on a machine with Docker.
python3 -m unittest discover -s tests -p test_tamarind.py -v
docker build -t zincsight-tamarind:test .
python3 tests/check_tamarind_container.py \
  --image zincsight-tamarind:test --output /tmp/zincsight-tamarind-results
```

The inherited preview includes reproducibility corrections that can change
historical predictions. Published accuracy benchmarks have not been re-established
for this revision. See [notebook regression](../mcp/NOTEBOOK_REGRESSION.md) and
[rotamer diagnosis](../mcp/ROTAMER_DIAGNOSIS.md). The adapter does not change
geometry/scoring formulas or trained weights. Predictions do not establish metal
identity, occupancy, or catalytic activity, and absence of a prediction does not
rule out binding.

Please cite [Mechtinger et al., Protein Science (2025)](https://doi.org/10.1002/pro.70350).
The existing repository license applies; see [LICENSE](../../LICENSE).

