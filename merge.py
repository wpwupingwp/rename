#!/usr/bin/python3

from os import mkdir, sep
from os.path import join as path_join
from os.path import basename
from glob import glob


folders = glob('*')
folders.remove('merge.py')
file_list = glob('*'+sep+'*')
file_list = [basename(i) for i in file_list]
choose = set(file_list)

mkdir('out')
for fasta in choose:
    to_merge = str()
    for folder in folders:
        try:
            with open(path_join(folder, fasta), 'r') as others:
                to_merge = ''.join([to_merge, others.read()])
        except:
            print('Miss', folder, fasta)
    with open(path_join('out', fasta), 'w') as out:
        out.write(to_merge)
