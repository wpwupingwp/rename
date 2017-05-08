#!/usr/bin/python3

from glob import glob
from os import mkdir
from os.path import join as join_path
from timeit import default_timer as timer
from Bio import SeqIO
# from sys import argv


def get_format():
    OLD = '1.gene|2.order|3.family|4.genus|5.species|6.accession_id|7.specimen'
    print(OLD)
    choice = input('Input number of field you want to use to divide:\n')
    return int(choice)


def divide(file_list, choice, out):
    for fasta in file_list:
        for record in SeqIO.parse(fasta, 'fasta'):
            item = record.id.split('|')[choice-1]
            # remove illegal characters
            item = item.replace('/', '_')
            handle = open(join_path(out, item+'.fasta'), 'a')
            SeqIO.write(record, handle, 'fasta')


def main():
    # gb = argv[1]
    file_list = glob('*.fasta')
    OUT = 'divide'
    mkdir(OUT)

    new_format = get_format()
    start = timer()
    divide(file_list, new_format, OUT)

    end = timer()
    print('''\nFinished with {0:.3f} s. You can find fasta file in the folder
    {1}.'''.format(end-start, OUT))


if __name__ == '__main__':
    main()
