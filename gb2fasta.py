#!/usr/bin/python3

from Bio import SeqIO
from os import mkdir
from os.path import join as join_path
from timeit import default_timer as timer
from sys import argv

from gene_rename import normalize


start = timer()
gb = argv[1]
# gb = input('Input file name:\n')
out = gb+'_out'
try:
    mkdir(out)
except:
    raise Exception('{0} already exists, please use another name.'.format(gb))
for record in SeqIO.parse(gb, 'gb'):
    """
    From Zhang guojin
    order statswith ales
    family startswith aceae except 8
    http://duocet.ibiodiversity.net/index.php?title=%E4%BA%92%E7%94%A8%E5%90%8D
    %E7%A7%B0&mobileaction=toggle_view_mobile"""
    # order|family|organims(genus|species)
    order_family = record.annotations['taxonomy']
    family_exception_raw = (
        'Umbelliferae,Palmae,Compositae,Cruciferae,Guttiferae,Leguminosae,'
        'Leguminosae,Papilionaceae,Labiatae,Gramineae')
    family_exception = family_exception_raw[0].split(',')
    order = ''
    family = ''
    for item in order_family:
        if item.endswith('ales'):
            order = item
            continue
        if item.endswith('aceae'):
            family = item
        elif item in family_exception:
            family = item
    organism = record.annotations['organism'].replace(' ', '_')
    genus, *species = organism.split('_')
    taxon = '{}|{}|{}|{}'.format(order, family, genus, '_'.join(species))

    accession = record.annotations['accessions'][0]
    try:
        specimen = record.features[0].qualifiers['specimen_voucher'
                                                 ][0].replace(' ', '_')
    except:
        specimen = ''
    for feature in record.features:
        if feature.type == 'gene' and 'gene' in feature.qualifiers:
            gene = feature.qualifiers['gene'][0].replace(' ', '_')
            gene = normalize(gene)[0]
            sequence = record.seq[feature.location.start:feature.location.end]
            sequence = str(sequence)
            with open(join_path(out, gene+'.fasta'), 'a') as handle:
                handle.write('>{}|{}|{}|{}\n{}\n'.format(
                    gene, taxon, accession, specimen, sequence))
end = timer()
print('''\nFinished with {0:.3f} s. You can find fasta file in the folder
{1}.'''.format(end-start, out))
