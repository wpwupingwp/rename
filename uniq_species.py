#!/usr/bin/python3

from Bio import SeqIO
from sys import argv
from timeit import default_timer as timer


def main():
    start = timer()
    raw = SeqIO.parse(argv[1], 'fasta')
    name_seq = dict()
    output = list()
    before = 0
    after = 0
    for record in raw:
        # gene|order|family|3.genus|4.species|accession_id|specimen
        name = record.id.split('|')[3:5]
        name = '.'.join(name)
        if name in name_seq:
            name_seq[name].append(record)
        else:
            name_seq[name] = [record, ]
        before += 1

    def get_longest(records):
        info = list()
        for record in records:
            n = record.seq.count('N') + record.seq.count('n')
            length = len(record.seq) - n
            info.append([length, record])
        longest = max(info, key=lambda x: x[0])
        return longest[1]

    print('Duplicated sequences:')
    log = open(argv[1]+'.log', 'w')
    for record in name_seq.values():
        after += 1
        longest = get_longest(record)
        if len(record) != 1:
            id_list = [i.id for i in record]
            log.write('Longest:\t{} in ({})\n'.format(
                longest.id, '\t'.join(id_list)))
    log.write('Before\t{}\tAfter\t{}\n'.format(before, after))
    with open(argv[1]+'.uniq', 'a') as output_file:
        SeqIO.write(output, output_file, 'fasta')
    end = timer()
    print('Total {} sequences in {} format.'.format(before, 'fasta'))
    print('{} sequences left in the file {}.uniq.'.format(after, argv[1]))
    print('Cost {:.3f} seconds.'.format(end-start))


if __name__ == '__main__':
    main()
