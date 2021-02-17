import sys, os, glob, re
import argparse
from Bio import SeqIO
from datetime import datetime

'''
usage example
 python3 scripts/pita_sequence_alignment_table.py -g gene_list -d pita_results_15-jan -r /home/lchapado/bio_programs/PITA/3UTRhuman.fa -m input_mirna -o aligment_result_minus_10_gene_list.tsv


'''

def check_arg (args=None) :
    '''
    The function is used for parsing the input parameters form the command line using the standard python package argparse.
    The package itself is handling  the validation and the return errors messages
    Input:
        args    # Contains the arguments from the command line
    Return:
        parser.parse_args()     # The variable contains the valid parameters
    '''
    parser = argparse.ArgumentParser(prog = 'summary_pita_results.py',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description= 'part of miRNA set of tools')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')

    #parser.add_argument('-g', '--genefile', required = False , help = 'file having the list of genes')

    parser.add_argument('-d', '--dir', required = True, help = 'directory with the PITA result files')

    parser.add_argument('-r', '--reference', required = True, help ='3UTR genome reference file')

    #parser.add_argument('-m', '--mirnas', required = True, help='directory with the miRNAs fasta files')

    parser.add_argument('-p', '--prefix', required = True, help='prefix file for output results')

    parser.add_argument('-e', '--score', required = True, help = 'filter by maximum energy score')

    return parser.parse_args()
