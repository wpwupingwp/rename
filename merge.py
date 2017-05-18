#!/usr/bin/python3

from os import mkdir, sep
from os.path import join as path_join
from os.path import basename
from glob import glob


folders = glob('*')
folders.remove('merge.py')
file_list = list()
for folder in folders:
    file_list.append([folder, glob(folder+sep+'*')])
choose = max(file_list, key=lambda i: len(i[1]))
folder = choose[0]
file_list = choose[1]
folders.remove(choose[0])

mkdir('out')
for fasta in file_list:
    with open(fasta, 'r') as head:
        to_merge = head.read()
    for folder in folders:
        try:
            with open(path_join(folder, basename(fasta)), 'r') as others:
                to_merge = ''.join([to_merge, others.read()])
        except:
            print('Miss', folder, fasta)
    with open(path_join('out', basename(fasta)), 'w') as out:
        out.write(to_merge)
