#!/usr/bin/env python3
"""Record implementation fingerprints independently of cached input validation."""
from pathlib import Path
import datetime,hashlib,json,sys
root=Path(sys.argv[1]).resolve()
paths=[p for folder in ['bin','modules','subworkflows','conf','vendor/bushman']
       for p in (root/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
paths += [root/'main.nf',root/'nextflow.config',root/'pipeline_identity.json',root/'vendor/bushman_r_sources/runtime_lock.json']
identity=json.loads((root/'pipeline_identity.json').read_text())
payload=dict(pipeline=identity['pipeline'],version=identity['version'],recorded_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    root=str(root),caller_implementation=sys.argv[2],site_window_bp=int(sys.argv[3]),
    abundance_method=sys.argv[4],annotation_mode=sys.argv[5],
    source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)})
Path('pipeline_implementation.json').write_text(json.dumps(payload,indent=2)+'\n')
