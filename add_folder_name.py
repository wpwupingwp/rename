#!/usr/bin/python3

from pathlib import Path

for folder in Path().glob('*'):
    if not folder.is_dir():
        continue
    name = folder.name
    for f in folder.glob('*'):
        f.rename(f.with_name(f'{name}-{f.name}'))
