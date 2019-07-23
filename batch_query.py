#!/usr/bin/python3

from Bio import Entrez
from time import sleep
import argparse
import os


def get_list(list_file):
    down_list = list()
    with open(list_file, 'r') as raw:
        for line in raw:
            down_list.append(line.strip())
    return down_list


def down(taxon_name, out_path):
    EMAIL = 'wpwupingwp@outlook.com'
    Entrez.email = EMAIL
    FILTER = '(plastid[filter] OR chloroplast[filter])'
    # set 'noexp' to fetch only this level
    query = '''"{}"[Organism] AND {} AND ("{}"[SLEN] : "{}"[SLEN])'''.format(
        taxon_name, FILTER, arg.min_len, arg.max_len)
    print(query)
    handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                        usehistory='y',))
    count = handle['Count']
    count = int(count)
    print('Totally {} records.'.format(count))
    if arg.skip:
        return
    file_name = os.path.join(out_path, taxon_name+'.gb')
    output_file = open(file_name, 'w')
    Retstart = 0
    while Retstart <= count:
        print(Retstart)
        try:
            genome_content = Entrez.efetch(db='nuccore',
                                           webenv=handle['WebEnv'],
                                           query_key=handle['QueryKey'],
                                           rettype='gb',
                                           retmode='text',
                                           retstart=Retstart,
                                           retmax=1000)
            output_file.write(genome_content.read())
        except KeyboardInterrupt:
            break
        except IOError:
            continue
        Retstart = Retstart + 1000
    print('{} finished.'.format(taxon_name))


def main():
    global arg
    arg = argparse.ArgumentParser()
    arg.add_argument('list_file',
                     help='Taxonomy name list, seperate by line')
    arg.add_argument('-min', dest='min_len', default=100000, type=int,
                     help='minium length')
    arg.add_argument('-max', dest='max_len', default=300000, type=int,
                     help='maximum length')
    arg.add_argument('-o', '-out', dest='out', default='Result',
                     help='output path')
    arg.add_argument('-skip', action='store_true',
                     help='only show records numbers')
    arg = arg.parse_args()
    if arg.skip is not True:
        os.mkdir(arg.out)
    for taxon in get_list(arg.list_file):
        down(taxon, arg.out)
        # Entrez limit
        sleep(0.5)


if __name__ == '__main__':
    main()
