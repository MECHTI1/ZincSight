"""Build the reviewable skills-only ZincSight plugin (Agent Plugins and Claude manifests) from an explicit file list."""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PLUGIN = ROOT / 'plugins/zincsight-preview'
FILES = (
    'plugin.json', '.claude-plugin/plugin.json', 'README.md', 'assets/icon.svg',
    'skills/zincsight-preview/SKILL.md',
    'skills/zincsight-preview/agents/openai.yaml',
    'skills/zincsight-preview/scripts/predict.py',
    'skills/zincsight-preview/scripts/requirements-lock.txt',
    'skills/zincsight-preview/assets/icon.svg',
)


def validate_listing(interface):
    """Check the fields this skills-only package uses against submission limits."""
    limits = {
        'displayName': 30, 'shortDescription': 30,
        'longDescription': 4000, 'developerName': 80,
    }
    for field, limit in limits.items():
        value = interface.get(field)
        if not isinstance(value, str) or not 0 < len(value.strip()) <= limit or len(value) > limit:
            raise ValueError(f'{field} must be nonempty and at most {limit} characters')
    prompts = interface['defaultPrompt']
    if not isinstance(prompts, list) or not 1 <= len(prompts) <= 3:
        raise ValueError('defaultPrompt must contain one to three starter prompts')
    if any(not isinstance(p, str) or not 0 < len(p.strip()) <= 128 or len(p) > 128 for p in prompts):
        raise ValueError('Starter prompts must be nonempty and at most 128 characters')
    if len(set(prompts)) != len(prompts):
        raise ValueError('Starter prompts must be unique')
    for field in ('composerIcon', 'logo'):
        name = interface[field]
        if not name.startswith('./') or name[2:] not in FILES:
            raise ValueError(f'{field} must reference a packaged relative path')
        path = PLUGIN / name[2:]
        # This package ships SVG artwork. Numeric dimensions must meet the
        # listing minimum even when the internal viewBox uses smaller units.
        if path.suffix != '.svg':
            raise ValueError(f'{field}: this builder currently validates SVG icons only')
        root = ET.parse(path).getroot()
        width, height = float(root.get('width', 0)), float(root.get('height', 0))
        if width != height or width < 48:
            raise ValueError(f'{field} must have square numeric dimensions of at least 48 pixels')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'dist')
    args = parser.parse_args()
    manifest = json.loads((PLUGIN / 'plugin.json').read_text())
    assert manifest['name'] == 'zincsight-preview'
    assert manifest['$schema'] == 'https://agent-plugins.org/schemas/1.0.0/plugin.schema.json'
    claude = json.loads((PLUGIN / '.claude-plugin/plugin.json').read_text())
    # Claude reads its own manifest; keep the shared identity fields in step.
    for field in ('name', 'version', 'description', 'author', 'homepage', 'repository', 'keywords'):
        assert claude.get(field) == manifest.get(field), f'.claude-plugin/plugin.json {field} differs from plugin.json'
    interface = manifest['extensions']['com.openai']['interface']
    validate_listing(interface)
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
