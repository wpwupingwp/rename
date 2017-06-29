#!/usr/bin/python3

from Bio import Entrez
from timeit import default_timer as timer
import argparse


def get_id_list(id_list_file):
    N = 1000
    with open(id_list_file) as raw:
        to_down = list()
        for index, line in enumerate(raw):
            to_down.append(line.strip())
            if index % N == 0:
                yield to_down
                to_down = list()


def down(id_list_file, email, output):
    Entrez.email = email
    output_file = open(output, 'w')
    for to_down in get_id_list(id_list_file):
        id_list = ', '.join(to_down)
        print('id list', id_list)
        handle = Entrez.read(Entrez.esearch(db='nuccore',
                                          term=id_list,
                                          usehistory='y',))
        print('handle', handle)
        genome_content = Entrez.efetch(db='nuccore',
                                       webenv=handle['WebEnv'],
                                       query_key=handle['QueryKey'],
                                       rettype='gb',
                                       retmode='text')
        output_file.write(genome_content.read())


def parse_args():
    arg = argparse.ArgumentParser()
    arg.add_argument('id_list', help='Accession ID list downloaded from NCBI')
    arg.add_argument('-email', default='wpwupingwp@outlook.com',
                     help='email address')
    arg.add_argument('-out', default='sequence.gb', help='output filename')
    return arg.parse_args()


def main():
    start = timer()
    arg = parse_args()
    print(arg)

    down(arg.id_list, arg.email, arg.out)
    end = timer()
    print('Cost {:.3f} seconds.'.format(end-start))


if __name__ == '__main__':
    main()
