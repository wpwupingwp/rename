#!/usr/bin/python3

import argparse
import re
from timeit import default_timer as timer
from subprocess import run
from os import mkdir, remove, cpu_count
from os.path import join as join_path
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from gene_rename import rename


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    arg.add_argument('gbfile', help='gb format file')
    arg.add_argument('-out',  help='output directory')
    arg.add_argument('-no_rename', action='store_true',
                     help='try to rename gene')
    arg.add_argument('-align', action='store_true',
                     help='use mafft to alignment')
    return arg.parse_args()


def safe(old):
    return re.sub(r'\W', '_', old)


def get_taxon(lineage, family_exception):
    superkingdom = ''
    kingdom = None
    order = ''
    family = ''
    for item in lineage:
        if item.endswith('ales'):
            order = item
        elif (item.endswith('aceae') or
              # elif (item.endswith('aceae') or item.endswith('idae') or
              item.endswith('viridae')):
            family = item
        elif item in family_exception:
            family = item
        elif item in ('Bacteria', 'Archaea', 'Viruses',
                      'Eukaryota', 'Viroids'):
            superkingdom = item
        elif item in ('Metazoa', 'Fungi', 'Viridiplantae'):
            kingdom = item
    if superkingdom == 'Eukaryota' and kingdom is not None:
        Type = kingdom
    else:
        Type = superkingdom
    return Type, order, family


def write_seq(name, sequence_id, feature, whole_seq, path, arg):
    """
    Write fasta file.
    """
    filename = join_path(path, name+'.fasta')
    sequence = feature.extract(whole_seq)

    with open(filename, 'a') as handle:
        handle.write(sequence_id+'\n')
        handle.write(str(sequence)+'\n')
    return filename


def get_feature_name(feature, arg):
    """
    Get feature name and collect genes for extract spacer.
    Only handle gene, product, misc_feature, misc_RNA.
    Return: [name, feature.type]
    """
    name = None
    misc_feature = None
    gene = None
    if feature.type == 'gene':
        if 'gene' in feature.qualifiers:
            gene = feature.qualifiers['gene'][0].replace(' ', '_')
            if not arg.no_rename:
                gene = rename(gene)[0]
            name = safe(gene)
        elif 'product' in feature.qualifiers:
            product = feature.qualifiers['product'][0].replace(' ', '_')
            name = safe(product)
        elif 'locus_tag' in feature.qualifiers:
            locus = feature.qualifiers['locus_tag'][0].replace(' ', '_')
            name = safe(locus)
        else:
            pass
    elif feature.type == 'misc_feature':
        misc_feature = None
        if 'product' in feature.qualifiers:
            misc_feature = feature.qualifiers['product'][0].replace(
                ' ', '_')
        elif 'note' in feature.qualifiers:
            misc_feature = feature.qualifiers['note'][0].replace(
                ' ', '_')
        if (misc_feature is not None) and ('intergenic_spacer' in misc_feature
                                           or 'IGS' in misc_feature):
            # 'IGS' in misc_feature) and len(misc_feature) < 100):
            name = safe(misc_feature)
            name = name.replace('intergenic_spacer_region',
                                'intergenic_spacer')
    elif feature.type == 'misc_RNA':
        if 'product' in feature.qualifiers:
            misc_feature = feature.qualifiers['product'][0].replace(
                ' ', '_')
        elif 'note' in feature.qualifiers:
            misc_feature = feature.qualifiers['note'][0].replace(
                ' ', '_')
        if misc_feature is not None:
            name = safe(misc_feature)
        else:
            return None, None
        # handle ITS
        if 'internal_transcribed_spacer' in name:
            name = 'ITS'
        # name = name.replace('internal_transcribed_spacer', 'ITS')
        # if 'ITS_1' in name:
        #     if 'ITS_2' in name:
        #         name = 'ITS'
        #     else:
        #         name = 'ITS_1'
        # elif 'ITS_2' in name:
        #     name = 'ITS_2'
    else:
        pass
        # print('Cannot handle feature:')
        # print(feature)
    return name, feature.type


def get_spacer(genes, arg):
    """
    List: [[name, SeqFeature],]
    """
    spacers = list()
    # sorted according to sequence starting postion
    genes.sort(key=lambda x: int(x[1].location.start))
    for n, present in enumerate(genes[1:], 1):
        before = genes[n-1]
        # use sort to handle complex location relationship of two fragments
        location = [before[1].location.start, before[1].location.end,
                    present[1].location.start, present[1].location.end]
        location.sort(key=lambda x: int(x))
        start, end = location[1:3]
        if before[1].location.strand == present[1].location.strand == -1:
            strand = -1
        else:
            strand = 1
        name = '_'.join([before[0], present[0]])
        spacer = SeqFeature(FeatureLocation(start, end), id=name,
                            type='spacer', strand=strand)
        spacers.append(spacer)
    return spacers


def divide(arg):
    """
    Given genbank file, return divided fasta files.
    From Zhang guojin
    order end with ales
    family end with aceae except 8
    http://duocet.ibiodiversity.net/index.php?title=%E4%BA%92%E7%94%A8%E5%90%8D
    %E7%A7%B0&mobileaction=toggle_view_mobile
    """
    # kingdom|order|family|organims(genus|species)
    family_exception_raw = (
        'Umbelliferae,Palmae,Compositae,Cruciferae,Guttiferae,Leguminosae,'
        'Leguminosae,Papilionaceae,Labiatae,Gramineae')
    family_exception = family_exception_raw[0].split(',')
    start = timer()
    groupby_gene = join_path(arg.out, '{}-groupby_gene'.format(arg.out))
    mkdir(groupby_gene)
    groupby_name = join_path(arg.out, '{}-groupby_name'.format(arg.out))
    mkdir(groupby_name)
    handle_raw = open(arg.gbfile+'.fasta', 'w')
    wrote_by_gene = set()
    wrote_by_name = set()

    for record in SeqIO.parse(arg.gbfile, 'gb'):
        # only accept gene, product, and spacer in misc_features.note
        lineage = record.annotations['taxonomy']
        # kingdom or superkingdom
        kingdom, order, family = get_taxon(lineage, family_exception)
        organism = record.annotations['organism'].replace(' ', '_')
        genus, *species = organism.split('_')
        taxon = '|'.join([kingdom, order, family, genus, '_'.join(species)])
        accession = record.annotations['accessions'][0]
        try:
            specimen = record.features[0].qualifiers['specimen_voucher'
                                                     ][0].replace(' ', '_')
        except (IndexError, KeyError):
            specimen = ''
        whole_seq = record.seq
        feature_name = list()
        genes = list()

        for feature in record.features:
            name, feature_type = get_feature_name(feature, arg)
            # skip unsupport feature
            if name is None:
                continue
            if len(name) > 100:
                print('Too long name: {}.'.format(name))
                name = name[:100] + '...'
            # skip abnormal annotation
            if len(feature) > 20000:
                print('Skip abnormal annotaion of {}!'.format(name))
                print('Accession: ', accession)
                continue
            if feature_type == 'gene':
                genes.append([name, feature])
            feature_name.append(name)
            sequence_id = '>' + '|'.join([name, taxon, accession, specimen])
            wrote = write_seq(name, sequence_id, feature, whole_seq,
                              groupby_gene, arg)
            wrote_by_gene.add(wrote)

        # extract spacer
        spacers = get_spacer(genes, arg)
        for spacer in spacers:
            sequence_id = '>' + '|'.join([spacer.id, taxon,
                                          accession, specimen])
            wrote = write_seq(spacer.id, sequence_id, spacer, whole_seq,
                              groupby_gene, arg)
            wrote_by_gene.add(wrote)
        # write to group_by name, i.e., one gb record one fasta
        if 'ITS' in feature_name:
            name_str = 'ITS'
        elif len(feature_name) >= 4:
            name_str = '{}-...-{}'.format(feature_name[0], feature_name[-1])
        elif len(feature_name) == 0:
            name_str = 'Unknown'
        else:
            name_str = '-'.join(feature_name)
        record.id = '|'.join([name_str, taxon, accession, specimen])
        record.description = ''
        filename = join_path(groupby_name, name_str+'.fasta')
        with open(filename, 'a') as out:
            SeqIO.write(record, out, 'fasta')
            wrote_by_name.add(filename)
        # write raw fasta
        SeqIO.write(record, handle_raw, 'fasta')

    end = timer()
    print('Divide done with {:.3f}s.'.format(end-start))
    return wrote_by_gene, wrote_by_name


def mafft(files):
    failed = list()
    # get available CPU cores
    cores = cpu_count()
    print('Start mafft ...')
    for fasta in files:
        print('Aligning {}'.format(fasta))
        out = fasta + '.aln'
        _ = ('mafft --thread {} --reorder --quiet --adjustdirection '
             ' {} > {}'.format(cores-1, fasta, out))
        m = run(_, shell=True)
        if m.returncode != 0:
            failed.append(out)
    print('Done with mafft.')
    return failed


def main():
    arg = parse_args()
    if arg.out is None:
        arg.out = arg.gbfile.replace('.gb', '')
    mkdir(arg.out)
    wrote_by_gene, wrote_by_name = divide(arg)
    if arg.align:
        failed = mafft(wrote_by_gene)
        for i in failed:
            print('Remove empty (failed alignment) file {}.'.format(i))
            remove(i)
    return


if __name__ == '__main__':
    main()
