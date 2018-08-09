##### 09.08.2018 #####
##### by Kevin Schneider #####
##### KevinSchneider@gmx.at #####
##### written for Python 3 #####
"""
--> script to filter fasta files <--
depending on the chosen flags/parameters, it filters:
     1) duplicates
     2) by the number of sequence records
     3) by the number of taxa
"""

from glob import glob
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
from os import path
import argparse

def get_unique(records):
     """ check if sequence is present more than once in the list of sequences """
     """ instead of comparing the sequence records explicitly, checksums are applied """
     checksums = set()
     for record in records:
          checksum = seguid(record.seq)
          if checksum in checksums:
               continue
          checksums.add(checksum)
          yield record
          
def check_rec_number(records, rec_l_limit, rec_u_limit):
     """ filter by the number of records in the file """
     rec_num = len(records)
     if rec_num > rec_l_limit and rec_num < rec_u_limit:
          return True

def check_taxon_number(records, taxa_limit):
     """ filter by the number of distinct taxa in the file """
     IDs = set()
     taxa_num = 0
     for record in records:
          if not record.id in IDs:
               IDs.add(record.id)
               taxa_num += 1
     if taxa_num > taxa_limit:
          return True

""" handling and parsing the input parameters """     
parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True, type=str, help='input files in fasta format (in quotes if wildcard [*] included)')
parser.add_argument('-d', action='store_true', help='remove duplicates')
parser.add_argument('-rl', nargs=1, type=int, help='lower limit for number of records')
parser.add_argument('-ru', nargs=1, type=int, help='upper limit for number of records')
parser.add_argument('-t', nargs=1, type=int, help='lower limit for number of taxa')
args = parser.parse_args()
files = args.i
remove_dup = args.d
rec_l_limit = args.rl[0]
rec_u_limit = args.ru[0]
taxa_limit = args.t[0]

""" loop through all the files with the specified file name (pattern) and open them in turn """  
for file in glob(files):
     with open(file, "rU") as f:
          """ create a list of sequence records in the file """
          records = list(SeqIO.parse(f, "fasta"))
          if remove_dup:
               new_records = get_unique(records)
          if not check_rec_number(records, rec_l_limit, rec_u_limit):
               continue
          if not check_taxon_number(records, taxa_limit):
               continue
          base = path.splitext(file)[0]
          ext = path.splitext(file)[1]
          new_file = base + ".filtered" + ext
          SeqIO.write(new_records, new_file, "fasta")
     
     
     
     

          
