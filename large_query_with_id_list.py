#!/usr/bin/python3

from Bio import Entrez


EMAIL = 'wpwupingwp@outlook.com'
Entrez.email = EMAIL
# set 'noexp' to fetch only this level
query = input('Query string:\n')
print(query)
handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                    usehistory='y',))
count = handle['Count']
Retstart = 1
output_file = open('sequence.gb', 'w')
while True:
    print(Retstart)
    genome_content = Entrez.efetch(db='nuccore',
                                   webenv=handle['WebEnv'],
                                   query_key=handle['QueryKey'],
                                   rettype='gb',
                                   retmode='text',
                                   retstart=Retstart,
                                   Retmax=100000)
    output_file.write(genome_content.read())
    Retstart = Retstart + 100000
