# Sharing ZincSight Preview

The simplest assistant package is one skill, running embedded PostgreSQL inside
the task. No hosted MCP server, paid database or predictor API key is needed.
The host's execution access and usage limits still apply. The original
[public Colab notebook](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/master/ZincSight.ipynb)
remains the browser fallback.

This is an experimental integration. The [diagnosis](ROTAMER_DIAGNOSIS.md)
documents known scientific-output differences. Keep the Preview name and explicit
invocation while scientific review remains incomplete.

**Package test passed on 1 November 2026:**
[run 41299312452](https://github.com/MECHTI1/ZincSight/actions/runs/41299312452)
built the archive, extracted it into an empty directory, installed dependencies
with a fresh cache, fetched public PDB 1KLS and returned its expected coordinating
residues with embedded PostgreSQL and the notebook defaults. An invalid file
request correctly produced no prediction. This validates package execution, not
a new biological benchmark or public-directory approval.

## Files already prepared

- `plugins/zincsight-preview/plugin.json`: portable Agent Plugins manifest,
  author information, discoverability keywords and three starter prompts.
- `plugins/zincsight-preview/skills/zincsight-preview/`: the tested skill,
  pinned launcher and dependency lock, exported from the installed preview.
- `.agents/plugins/marketplace.json`: a repository catalog for supported clients.
- `tools/build_plugin.py`: an explicit-file-list package builder, excluding all
  private inputs, outputs, credentials, environments and historical archives.
- `tools/check_plugin_install.py`: extraction and real launcher check, with a
  fresh dependency cache by default and a public 1KLS control.
- The **Plugin package** GitHub Actions workflow builds the archive, runs that
  check, and retains the package, checksum and result as one workflow artifact.

## Try through a repository marketplace

For clients supporting GitHub marketplaces, add this preview branch:

```sh
codex plugin marketplace add MECHTI1/ZincSight --ref feature/mcp-integration
```

Open the supported client's Plugins Directory, select the ZincSight Preview
marketplace, and install ZincSight Preview. Marketplace import availability
depends on the client and workspace. A branch import is a preview test, not a
public directory listing. To use a reviewed immutable snapshot, replace the
branch argument with its commit SHA.

Ask:

> Use ZincSight Preview on PDB 1KLS. Show the coordinating residues, predicted
> zinc position and probability. Save the full JSON.

The launcher installs its own pinned environment on first use. It needs Python
3.12, Node.js 22+, npm, a C compiler, PyMOL system libraries and allowed downloads.
If those are blocked, use Colab; changing the database does not grant network or
runtime permissions.

## Prepare the public directory submission

Build locally with `python tools/build_plugin.py`, or obtain the tested archive
from the **zincsight-preview-package** artifact on a successful **Plugin package**
workflow run. The inner file `zincsight-preview-0.1.0-preview.1.zip` is the plugin
package; the outer GitHub artifact download is only a container for it.

In OpenAI's plugin submission portal, choose the skills-only route, upload the
inner plugin ZIP under the verified publisher identity, review scan findings,
and submit the draft when its scientific scope and required declarations have
been reviewed. Publishing follows approval. The account owner must supply or
complete identity verification if it is missing; this repository cannot do that.

This package intentionally contains no MCP connection, lifecycle hooks or
account-specific app references. Skills-only submissions do not require a
permanent HTTPS endpoint. Do not add the local stdio MCP server to this cloud
package: that has different distribution requirements. A future hosted MCP
integration needs its own supported submission path.

Installing a personal skill, sharing it within a workspace, committing this
package to GitHub, and being listed in the public directory are distinct steps.
None guarantees that every ChatGPT user will be shown or automatically use it.
The package remains explicit-use while labeled Preview.

## Text prepared for a preview invitation

> Try ZincSight Preview for predicting potential zinc-binding sites directly
> from protein structures. Provide a PDB ID, AlphaFold/UniProt identifier, MGYP
> ID or PDB/mmCIF file, and inspect coordinating residues, zinc positions and
> probability scores. The assistant integration runs on demand with an embedded
> database and requires no separate ZincSight API key. This is an experimental
> integration; results can differ from the original notebook. The original
> public Colab workflow is also available.

Use the repository or Colab URL when inviting testers. Add a public directory
URL only after one exists. This text is a draft; no announcement has been sent.
Do not apply the published ZincSight accuracy benchmarks to this revised engine
without scientific validation. No software license is selected by this package.

## Official distribution references

Checked 1 November 2026:

- [Plugin packaging and marketplace format](https://developers.openai.com/plugins/build/plugins)
- [Upload, review and publication](https://developers.openai.com/plugins/deploy/submission)
- [Skills-only submission guidance](https://developers.openai.com/plugins/guides/submit-claude-plugin)

## Maintaining the skill export

Keep the exported `SKILL.md`, `agents/openai.yaml`, launcher, dependency lock and
icon identical to the tested installed skill. For changes, first update and test
that skill, export only those five files, bump the plugin version, and rerun the
package check. The engine source remains pinned by revision and archive hash;
editing other branch files does not silently change what the plugin executes.
