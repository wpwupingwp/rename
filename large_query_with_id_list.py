#!/usr/bin/python3

from Bio import Entrez
from timeit import default_timer as timer
import argparse


def get_id_list(id_list_file, batch_size, redo):
    id_list = list()
    with open(id_list_file) as raw:
        for line in raw:
            id_list.append(line.strip())
    if redo is not None:
        try:
            location = id_list.index(redo)
            id_list = id_list[location:]
            print('Skipped {} records.'.format(location))
        except:
            raise ValueError('Wrong redo value! It looks like\n{}'.format(
                id_list[0]))
    for n in range(0, len(id_list), batch_size):
        # list slice return empty list if slice out of len(list)
        yield n, id_list[n:(n+batch_size)]


def down_wrapper(id_list_file, batch_size, email, output, redo):
    Entrez.email = email
    out = open(output, 'w')
    tried = 0
    last_one = ''
    for n, to_down in get_id_list(id_list_file, batch_size, redo):
        while True:
            if down(n, to_down, out):
                break
            elif tried == 3:
                print('Too much failure. Please check your network.')
                print('Retry ...')
            elif tried == 10:
                print('Seems not work. Quit now!')
                print('If you want to continue, add option'
                      ' "-redo {}" in old command'.format(last_one))
                break
            else:
                print('Retry ...')
            tried += 1
        last_one = to_down[-1]
    out.close()


def down(n, to_down, output_file):
    try:
        print('{}: {} ... {}'.format(n, to_down[0], to_down[-1]), end='\t')
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
        got = len(to_down)
        print('{} records got.'.format(got))
        return True
    except KeyboardInterrupt:
        raise
    except:
        return False


def parse_args():
    arg = argparse.ArgumentParser()
    arg.add_argument('id_list', help='Accession ID list downloaded from NCBI')
    arg.add_argument('-e', '--email', default='wpwupingwp@outlook.com',
                     help='email address')
    arg.add_argument('-b', '--batch_size', type=int, default=500,
                     help='Number of records to download once')
    arg.add_argument('-redo', help='id which restart from')
    arg.add_argument('-out', default='sequence.gb', help='output filename')
    return arg.parse_args()


def main():
    start = timer()
    global arg
    arg = parse_args()

    down_wrapper(arg.id_list, arg.batch_size, arg.email, arg.out, arg.redo)
    end = timer()
    print('Cost {:.3f} seconds.'.format(end-start))


if __name__ == '__main__':
    main()
