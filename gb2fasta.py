#!/usr/bin/python3

import re
from Bio import SeqIO
from os import mkdir
from os.path import join as join_path
from sys import argv
from timeit import default_timer as timer

from gene_rename import normalize


start = timer()

groupby_gene = '{}-groupby_gene'.format(argv[1].replace('.gb', ''))
mkdir(groupby_gene)
groupby_name = '{}-groupby_name'.format(argv[1].replace('.gb', ''))
mkdir(groupby_name)
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


def safe(old):
    return re.sub(r'\W', '_', old)


def get_seq(feature, whole_sequence, expand=False, expand_n=100):
    """
    Given location and whole sequence, return fragment.
    """
    if expand:
        # in case of negative start
        start = max(0, feature.location.start-expand_n)
        end = min(len(whole_sequence), feature.location.end+expand_n)
    else:
        start = feature.location.start
        end = feature.location.end
    return whole_sequence[start:end]


def extract(feature, name, whole_seq):
    filename = join_path(groupby_gene, name+'.fasta')
    sequence = get_seq(feature, whole_seq, expand=False)
    with open(filename, 'a') as handle:
        handle.write('>{}|{}|{}|{}\n{}\n'.format(
            name, taxon, accession, specimen, sequence))
    filename2 = join_path(groupby_gene, 'expand.{}.fasta'.format(name))
    sequence = get_seq(feature, whole_seq, expand=True, expand_n=100)
    with open(filename2, 'a') as handle:
        handle.write('>{}|{}|{}|{}\n{}\n'.format(
            name, taxon, accession, specimen, sequence))


for record in SeqIO.parse(argv[1], 'gb'):
    # only accept gene, product, and spacer in misc_features.note
    order_family = record.annotations['taxonomy']
    order, family = get_taxon(order_family)
    organism = record.annotations['organism'].replace(' ', '_')
    genus, *species = organism.split('_')
    taxon = '{}|{}|{}|{}'.format(order, family, genus, '_'.join(species))
    accession = record.annotations['accessions'][0]
    try:
        specimen = record.features[0].qualifiers['specimen_voucher'
                                                 ][0].replace(' ', '_')
    except (IndexError, KeyError):
        specimen = ''
    whole_seq = record.seq
    feature_name = list()

    for feature in record.features:
        gene = ''
        prodcut = ''
        if feature.type == 'gene':
            if 'gene' in feature.qualifiers:
                gene = feature.qualifiers['gene'][0].replace(' ', '_')
                gene = normalize(gene)[0]
                name = safe(gene)
            elif 'product' in feature.qualifiers:
                product = feature.qualifiers['product'][0].replace(' ', '_')
                name = safe(product)
        elif feature.type == 'misc_feature' and 'note' in feature.qualifiers:
            misc_feature = feature.qualifiers['note'][0].replace(' ', '_')
            if (('intergenic_spacer' in misc_feature or 'IGS' in misc_feature)
                    and len(misc_feature) < 100):
                name = safe(misc_feature)
                name = name.replace('intergenic_spacer_region',
                                    'intergenic_spacer')
            else:
                print(misc_feature)
                continue
        else:
            continue
        extract(feature, name, whole_seq)
        feature_name.append(name)
    # if len(genes) > 4:
    #     name_str = '{}-{}-{}genes-{}-{}'.format(*genes[:2], len(genes)-4,
    #                                             *genes[-2:])
    # elif 0 < len(genes) <= 4:
    #     name_str = '-'.join(genes)
    # elif len(products) > 4:
    #     name_str = '{}-{}-{}products-{}-{}'.format(
    #         *products[:2], len(products)-4, *products[-2:])
    # elif 0 < len(products) <= 4:
    #     name_str = '-'.join(products)
    # elif len(misc_features) > 4:
    #     name_str = '{}-{}-{}misc_features-{}-{}'.format(
    #         *misc_features[:2], len(misc_features)-4, *misc_features[-2:])
    # elif 0 < len(misc_features) <= 4:
    #     name_str = '-'.join(misc_features)
    # else:
    #     name_str = 'Unknown'
    if len(feature_name) >= 4:
        name_str = '{}-{}features-{}'.format(
            feature_name[0], len(feature_name)-2, feature_name[-1])
    elif len(feature_name) != 0:
        name_str = '-'.join(feature_name)
    else:
        name_str = 'Unknown'

    record.id = '|'.join([name_str, taxon, accession, specimen])
    record.description = ''
    SeqIO.write(record, handle_raw, 'fasta')
    with open(join_path(groupby_name, name_str+'.fasta'), 'a') as handle_name:
        SeqIO.write(record, handle_name, 'fasta')


end = timer()
print('Done with {:.3f}s.'.format(end-start))
