#!/usr/bin/python3

from Bio import Entrez
from timeit import default_timer as timer
import argparse


def get_id_list(id_list_file):
    N = 500
    id_list = list()
    with open(id_list_file) as raw:
        for line in raw:
            id_list.append(line.strip())
    for n in range(0, len(id_list), N):
        # list slice return empty list if slice out of len(list)
        yield id_list[n:(n+N)]


def down(id_list_file, email, output):
    Entrez.email = email
    output_file = open(output, 'w')
    for to_down in get_id_list(id_list_file):
        print('{} ... {}'.format(to_down[0], to_down[-1]))
        id_list = ','.join(to_down)
        handle = Entrez.read(Entrez.esearch(db='nuccore',
                                            term=id_list,
                                            usehistory='y',))
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

    down(arg.id_list, arg.email, arg.out)
    end = timer()
    print('Cost {:.3f} seconds.'.format(end-start))


if __name__ == '__main__':
    main()
