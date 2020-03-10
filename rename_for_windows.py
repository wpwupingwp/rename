#!/usr/bin/python3

from pathlib import Path

d_suffix = '.delete_because_duplicate'
files = list(Path('.').glob('*'))
files_lower = dict()
for f in files:
    if f.suffix == d_suffix:
        continue
    lower = str(f).lower()
    if lower in files_lower:
        print(f, files_lower[lower])
        # pathlib do not support write append
        with open(files_lower[lower], 'a') as out:
            out.write(f.read_text())
        f.rename(f.with_suffix(d_suffix))
    else:
        files_lower[lower] = f
