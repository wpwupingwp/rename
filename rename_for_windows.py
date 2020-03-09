#!/usr/bin/python3

from pathlib import Path

files = list(Path('.').glob('*'))
files_lower = dict()
for f in files:
    lower = str(f).lower()
    if lower in files:
        # pathlib do not support write append
        with open(files_lower[lower], 'a') as out:
            out.write(f.read_text())
        f.rename(f.with_suffix('.delete'))
    files_lower[lower] = f
