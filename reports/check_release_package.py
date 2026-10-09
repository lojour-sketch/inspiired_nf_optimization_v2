#!/usr/bin/env python3
"""Check the portable source/report tree without executing scientific analysis."""
import hashlib
import json
import subprocess
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
files = sorted(p for p in ROOT.rglob('*') if '.git' not in p.relative_to(ROOT).parts and p.is_file())
failures = []
checks = {}

def check(label, ok, detail=''):
    checks[label] = {'passed': bool(ok), 'detail': detail}
    if not ok:
        failures.append(label + ': ' + detail)

for path in files:
    relative = str(path.relative_to(ROOT))
    check('regular_file:' + relative, not path.is_symlink(), 'No external symlinks in the code tree')
    check('size:' + relative, path.stat().st_size < 50 * 1024 * 1024, 'Code/report files must be below 50 MiB')
    forbidden = (path.suffix.lower() in ('.sif', '.2bit', '.rdata', '.psl', '.o', '.so') or
                 any(part in ('work', '.nextflow', '__pycache__') for part in path.relative_to(ROOT).parts) or
                 '.fastq' in path.name.lower() or '.fq.' in path.name.lower())
    check('source_only:' + relative, not forbidden, 'No raw reads, full scientific archives, runtime or work directory')
    if path.suffix == '.json':
        json.loads(path.read_text())
    if path.suffix == '.py' and 'vendor' not in path.relative_to(ROOT).parts:
        compile(path.read_text(), str(path), 'exec')
    if path.suffix == '.sh':
        result = subprocess.run(['bash', '-n', str(path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        check('shell_syntax:' + relative, result.returncode == 0, result.stderr.decode())

provenance = json.loads((ROOT / 'vendor/bushman/PROVENANCE.json').read_text())
for record in provenance['files']:
    path = ROOT / 'vendor/bushman' / record['component'] / record['file']
    check('upstream_hash:' + str(path.relative_to(ROOT)),
          hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256'], 'Frozen upstream SHA256')

class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.urls = []

    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name in ('href', 'src') and value:
                self.urls.append(value)

local_links = 0
for path in files:
    if path.suffix != '.html':
        continue
    parser = Links()
    parser.feed(path.read_text())
    for url in parser.urls:
        parsed = urlsplit(url)
        if parsed.scheme or parsed.netloc or not parsed.path:
            continue
        local_links += 1
        target = (path.parent / unquote(parsed.path)).resolve()
        check('html_link:' + str(path.relative_to(ROOT)) + ':' + url,
              target.exists() and ROOT in target.parents, 'Portable local HTML link')

summary = {'status': 'passed' if not failures else 'failed', 'files_checked': len(files),
           'local_html_links_checked': local_links,
           'upstream_files_verified': len(provenance['files']),
           'scope': 'Package contents, Python/shell syntax, JSON, frozen source hashes and portable report links; scientific tests recorded separately',
           'failures': failures}
(ROOT / 'reports/evidence/release_package_check.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
if failures:
    raise SystemExit(1)
