#!/usr/bin/env python3
import argparse
import json
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from prepare_fasta import sha256

p = argparse.ArgumentParser()
p.add_argument('--compiler', required=True)
p.add_argument('--command', required=True)
p.add_argument('sources', nargs='+')
a = p.parse_args()
info = dict(created_at=datetime.now(timezone.utc).isoformat(), command=a.command,
    compiler=subprocess.check_output(shlex.split(a.compiler) + ['--version'], text=True),
    binary_sha256=sha256('align_bench'),
    sources={s: sha256(s) for s in a.sources})
Path('build.json').write_text(json.dumps(info, indent=2) + '\n')
