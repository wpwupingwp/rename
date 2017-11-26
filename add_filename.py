#!/usr/bin/python3

from sys import argv

with open(argv[1], 'r') as old, open(argv[1]+'.new', 'w') as new:
    for line in old:
        if line.startswith('>'):
            new.write(line.replace('>', '>{}-'.format(argv[1])))
        else:
            new.write(line)
