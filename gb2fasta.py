#!/usr/bin/python3

import re
from Bio import SeqIO
from os import mkdir
from os.path import join as join_path
from timeit import default_timer as timer
from sys import argv

from gene_rename import normalize


def safe(old):
    return re.sub(r'\W', '_', old)


start = timer()
out = '{}-out'.format(argv[1].replace('.gb', ''))
mkdir(out)
handle_raw = open(argv[1]+'.fasta', 'w')


family_exception_raw = (
    'Umbelliferae,Palmae,Compositae,Cruciferae,Guttiferae,Leguminosae,'
    'Leguminosae,Papilionaceae,Labiatae,Gramineae')
family_exception = family_exception_raw[0].split(',')


def get_taxon(taxonomy):
    """
    From Zhang guojin
    order statswith ales
    family startswith aceae except 8
    http://duocet.ibiodiversity.net/index.php?title=%E4%BA%92%E7%94%A8%E5%90%8D
    %E7%A7%B0&mobileaction=toggle_view_mobile"""
    # order|family|organims(genus|species)
    order = ''
    family = ''
    for item in order_family:
        if item.endswith('ales'):
            order = item
        elif item.endswith('aceae'):
            family = item
        elif item in family_exception:
            family = item
    return order, family


for record in SeqIO.parse(argv[1], 'gb'):
    print(record.name)
    order_family = record.annotations['taxonomy']
    order, family = get_taxon(order_family)
    organism = record.annotations['organism'].replace(' ', '_')
    genus, *species = organism.split('_')
    taxon = '{}|{}|{}|{}'.format(order, family, genus, '_'.join(species))
    accession = record.annotations['accessions'][0]
    try:
        specimen = record.features[0].qualifiers['specimen_voucher'
                                                 ][0].replace(' ', '_')
    except KeyError:
        specimen = ''
    seq = record.seq
    for feature in record.features:
        if feature.type == 'gene' and 'gene' in feature.qualifiers:
            gene = feature.qualifiers['gene'][0].replace(' ', '_')
            gene = normalize(gene)[0]
            gene = safe(gene)
            try:
                sequence = feature.extract(seq)
            except ValueError:
                sequence = ''
            sequence = str(sequence)
            with open(join_path(out, gene+'.fasta'), 'a') as handle:
                handle.write('>{}|{}|{}|{}\n{}\n'.format(
                    gene, taxon, accession, specimen, sequence))
    record.description = ''
    record.id = '|'.join(['raw', taxon, accession, specimen])
    SeqIO.write(record, handle_raw, 'fasta')

end = timer()
print('''\nFinished with {0:.3f} s. You can find fasta file in the folder
{1}.'''.format(end-start, out))
