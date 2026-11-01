"""Build the reviewable skills-only ZincSight plugin from an explicit file list."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PLUGIN = ROOT / 'plugins/zincsight-preview'
FILES = (
    'plugin.json', 'README.md', 'assets/icon.svg',
    'skills/zincsight-preview/SKILL.md',
    'skills/zincsight-preview/agents/openai.yaml',
    'skills/zincsight-preview/scripts/predict.py',
    'skills/zincsight-preview/scripts/requirements-lock.txt',
    'skills/zincsight-preview/assets/icon.svg',
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'dist')
    args = parser.parse_args()
    manifest = json.loads((PLUGIN / 'plugin.json').read_text())
    assert manifest['name'] == 'zincsight-preview'
    assert manifest['$schema'] == 'https://agent-plugins.org/schemas/1.0.0/plugin.schema.json'
    interface = manifest['extensions']['com.openai']['interface']
    assert all(0 < len(prompt) <= 128 for prompt in interface['defaultPrompt'])
    assert len(interface['defaultPrompt']) <= 3
    for field in ('composerIcon', 'logo'):
        assert interface[field].startswith('./')
        assert interface[field][2:] in FILES
    actual = {p.relative_to(PLUGIN).as_posix() for p in PLUGIN.rglob('*') if p.is_file()}
    assert actual == set(FILES), f'Unexpected or missing package files: {actual ^ set(FILES)}'
    for name in FILES:
        path = PLUGIN / name
        assert not path.is_symlink() and path.stat().st_size < 1_000_000, name
    compile((PLUGIN / 'skills/zincsight-preview/scripts/predict.py').read_text(), 'predict.py', 'exec')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / f"{manifest['name']}-{manifest['version']}.zip"
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(FILES):
            entry = zipfile.ZipInfo(name, date_time=(2026, 1, 1, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_DEFLATED
            entry.external_attr = 0o100644 << 16
            archive.writestr(entry, (PLUGIN / name).read_bytes())
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    output.with_suffix('.zip.sha256').write_text(f'{digest}  {output.name}\n')
    print(json.dumps({'archive': str(output), 'sha256': digest, 'files': len(FILES)}))


if __name__ == '__main__':
    main()
