#!/usr/bin/python3

import json
from sys import argv

with open('./genus_dict.json', 'r') as _:
    name_dict = json.load(_)

old = open(argv[1], 'r')
new = open(argv[1]+'.rename', 'w')
for line in old:
    if line[0] == '>':
        (gene, kingdom, order, family, genus, species, accession,
         *_) = line.strip().split('|')
        species_clean = species.split('_')[0]
        if family == '':
            try:
                kingdom, order, family, genus = name_dict[
                    genus+'|'+species_clean]
                new.write('|'.join([gene, kingdom, order, family, genus,
                                    species, accession, *_])+'\n')
            except KeyError:
                new.write(line)
        else:
            new.write(line)
    else:
        new.write(line)
