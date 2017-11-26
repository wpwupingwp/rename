# Name format

Most of programs below use sequence id format like this:

>  1.gene|2.order|3.family|4.genus|5.species|6.accession_id|7.specimen_voucher

# gene_rename.py

Covert genbank format to fasta format with reformated id.

Normalize gene name by these rules:

- rRNA gene should start with "rrn"
- tRNA gene should look like "trnXnnn", and "X" is the letter of amino acids,
  "nnn" is 3 letters of codon. Note that it is tran-reversed.
-  for other gene names, characters not being alphabet will be removed, and
   suffix like number or letter of subunit will be capitalized.
- name which could not be recognized will do nothing

## Usage

This is a python3 function. It need *re* and *Biopython* to run.

To use it, follow this example:

> from gene_rename import normalize

> new_name, name_type = normalize(old_name)

Mention that all input and output is string. For *name_type*, it will be:

* bad_name
* suspicious_name
* tRNA
* rRNA
* normal

"bad_name" means it cannot be recognized and will be same as input.
"suspicious_name" means it its too long (longer than 15 characters).

# gb2fasta.py

## Usage

> python3 gb2fasta.py gb_file-name

Find output files in "file_name_out".

## Example

> python3 gb2fasta.py sequence.gb

# fasta_rename.py

## Name format

Sequence id follows this format:

>  1.gene|2.order|3.family|4.genus|5.species|6.accession_id|7.specimen_voucher

## Usage

- Put this program in same folder with fasta files you want to change id. 
- Double click to run. 
-  Input number of fields you want in the order you want. Notice that every field could only be used once.
-  Find output files in "renamed" folder.

# fasta_rename_advance.py

## Name format

*This program can be used in any name format.*

## Usage

- Put this program in same folder with fasta files you want to change id. 
- Double click to run. 
-  Input number of fields you want in the order you want and seperators you want to use. Notice that every
field could only be used once.
-  Find output files in "renamed" folder.

# group_by.py

Divide sequences in input fasta file into different file according to
seperator you choose. Make sure they have same format in sequence id.

If the field as seperator in sequence id does not exist, it will be put into
FAILED.fasta.

## Usage

- Run
> python3 group_by.py input_file
- Choose the field you want to use to be seperator
- If you know which field you want the you can run like this:
> python3 group_by.py input_file -c n

The n is the number of field.
If you do not set "-c", it will hint you to choose which field of sequence id
to be used as separator.

## Example

> python3 group_by.py cgl.fasta -c 4

# large_query.py

Input query string same as in NCBI Genbank and download large data.

## Usage

> python3 large_query.py

Then input query string.

# large_query_with_id_list.py

According to accession id list you downloaded from NCBI, download gb records.
Note that this program will not check output. So you have to verify the data
by self.

Usage:

> python3 large_query_with_id_list.py id_list -redo accession_number

The id_list is accession list file you downloaded before. If the download
process failed, program will quit and you can use "-redo" to continue from
given accession number.

# uniq.py

Only left one record per species by remove shorter sequence (consider 'N').

## Usage

> python3 uniq-species.py input.fasta -c choice.

See log file for detail output.

If you do not set "-c", it will hint you to choose which field of sequence id
was used to divide sequence.

## Example

> python3 uniq.py whole.fasta -c "4 5"

# replace.py

Replace BOP with info given in "info.csv".

## Usage

> python3 replace.py

Make sure you have info.csv in the same folder with replace.py and fasta files
end with ".fasta".

Only handle ">BOP000000".

# add_filename.py

Add filename into the head of sequence id.

## Usage

> python3 add_filename.py fastafile

# Requirement

1. [python3](https://www.python.org/downloads/)

Be sure to install python3 rather than python 2.7. Besides, to use subprocess.run, you would better install python **3.5** or above

2. [biopython](http://biopython.org/wiki/Download)
 
Be sure to download same version of python3
