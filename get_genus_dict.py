#!/usr/bin/python3

import json

result = dict()
with open('./taxonomy.csv', 'r') as _:
    for line in _:
        Type, order, family, genus, *species = line.strip().split(',')
        # ignore subspecies
        species = species[0]
        name = '{}|{}'.format(genus, species)
        result[name] = [Type, order, family, genus]
with open('genus_dict.json', 'w') as out:
    json.dump(result, out)
