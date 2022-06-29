#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Filipe Z. Dezordi; Yago Dias"
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Filipe Z. Dezordi; Yago Dias"
__email__ = "zimmer.filipe@gmail.com; yag.dias@gmail.com"
__date__ = "2021/03/15"
__username__ = "dezordi; yagdias"

# libraries
from os import execl
from Bio import Entrez
from Bio import SeqIO
import argparse
import re
import csv


# arguments
parser = argparse.ArgumentParser(
    description='This script retrieve a formatted table for using on EEfinder', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument(
    "-in", "--input", help="File with access IDs",  required=True)
parser.add_argument(
    "-em", "--email", help="Email for using on Entrez", required=True)
parser.add_argument("-key", "--APIkey",
                    help="Key for not limiting number attemps", required=True)
args = parser.parse_args()
sequences = args.input
email = args.email
APIkey = args.APIkey


def fasta_ids(input_file):
    ids_list = []
    for record in SeqIO.parse(input_file, 'fasta'):
        ids_list.append(record.id)
    return ids_list


def efetch_function(data):
    Entrez.email = email
    Entrez.api_key = APIkey
    handle = Entrez.efetch(db='protein', id=data,
                           retmode='xml', rettype='gb')
    return handle.read()


if __name__ == '__main__':
    ids_list = fasta_ids(sequences)
    row_list = []

    with open('output_file', 'w', newline='') as output:
        writer = csv.writer(output, delimiter=',')
        header = ["Acession", "Species", 'Genus',
                  'Family', 'Molecule_type', 'Protein', 'Host']
        writer.writerow(header)
        for acession in ids_list:
            xml = efetch_function(acession)
            xml = xml.decode("utf-8")
            taxid = str(re.findall(r'taxon:.*', xml))
            taxid = re.sub(r'\[.*:', '', taxid)
            taxid = re.sub(r'<.*\]', '', taxid)
            handle_tax_data = Entrez.efetch(
                db="Taxonomy", id=taxid, retmode="xml")
            tax_info = Entrez.read(handle_tax_data)
            species = genus = family = ''
            try:
                for key in tax_info[0]["LineageEx"]:
                    if 'family' in key.values():
                        family = key['ScientificName']
                    if 'genus' in key.values():
                        genus = key['ScientificName']
                    if 'species' in key.values():
                        species = key['ScientificName']
                print(
                    f"done organism taxonomy for {acession}:{taxid}", end='\n')
            except:
                print(
                    f'error in organism taxonomy step for {acession}:{taxid}',  end='\n')
            molecule_type = 'DNA'
            protein = str(re.findall(r'<GBSeq_definition>.*', xml))
            protein = re.sub(r"\['<GBSeq_definition>", '', protein)
            protein = re.sub(r'<.*\]', '', protein)
            host = 'None'
            if family == '':
                family = 'Unclassified'
            if genus == '':
                genus = 'Unclassified'
            if species == '':
                species = 'Unclassified'
            row_list.append([acession, species, genus, family,
                            molecule_type, protein, host])
        writer.writerows(row_list)
