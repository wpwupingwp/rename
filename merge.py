#!/usr/bin/python3

from os import mkdir, walk
from os.path import join as path_join


# os.walk [this, subfolder, files] -R
file_list = list((walk('./')))
folder_list = file_list[0][1:][0]
example = folder_list[0]
folder_list = folder_list[1:]
file_list = file_list[1][2:][0]
mkdir('merge')
for fasta in file_list:
    with open(path_join(example, fasta), 'r') as head:
        to_merge = head.read()
    for folder in folder_list:
        try:
            with open(path_join(folder, fasta), 'r') as others:
                to_merge = ''.join([to_merge, others.read()])
        except:
            continue
    with open(path_join('merge', fasta), 'w') as out:
        out.write(to_merge)
