import sys, os, glob, re
import pandas as pd
import argparse
from Bio import SeqIO
from datetime import datetime
from utils.generic_functions import *

'''

usage example
 python3 ~/tesis_code/miRNAs_utils/summary_pita_results.py  -d pita_results_15-jan -r /home/lchapado/bio_programs/PITA/3UTRhuman.fa -p aligment -e -10

ath_miR165a_5p_pita_results.tab
'''
pita_result_estension = '_results.tab'
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


def fetch_target_genes_data(in_dir, filter_name, min_energy):
    '''
    Description:
        The function get the target genes form the PITA output files and reduce the
        results whith the ones that have less energy value as indicated in min_energy
        and it stores in
    Input:
        in_dir   # output PITA directory
        filter_name     # target to match the PITA output files
        min_energy      # minimun value to filter the target genes
    Return:
        target_genes_dict
    UTR	microRNA	Start	End	Seed	Loop	dGduplex	dG5	dG3	dG0	dG1	dGopen	ddG
    ENST00000310015|SP3	ath_miR165a-5p	6872	6864	8:0:1	0	-21.1	-7.1	-14	-29.49	-27.61	-1.87	-19.22
    ENST00000320285|AGPAT4	ath_miR165a-5p	4781	4773	8:0:1	0	-24.52	-9.3	-15.22	-41.37	-35.52	-5.85	-18.66
    '''
    target_genes_dict = {}
    gene_start_end_dict = {}
    filter_files = os.path.join(in_dir,'*' + filter_name)
    files_to_process =glob.glob(filter_files)
    for file in files_to_process:
        panda_df = pd.read_csv(file, sep ='\t', header=0)
        selected_df = panda_df[panda_df['ddG'] <= min_energy]
        selected_df[['UTR','Gene']] = selected_df['UTR'].str.split('|',expand = True)
        selected_df['Start_End'] = selected_df['Start'].apply(str).str.cat(selected_df['End'].apply(str),sep = '_')
        for index, row in selected_df.iterrows():
            if not row['Gene'] in target_genes_dict:
                target_genes_dict[row['Gene']] = {}
                target_genes_dict[row['Gene']][row['UTR']] ={}
                target_genes_dict[row['Gene']][row['UTR']][row['Start_End']] ={}
                target_genes_dict[row['Gene']][row['UTR']][row['Start_End']]['microRNA'] = row['microRNA']
                target_genes_dict[row['Gene']][row['UTR']][row['Start_End']]['Seed'] = row['Seed']
                target_genes_dict[row['Gene']][row['UTR']][row['Start_End']]['ddG'] = row['ddG']
                gene_start_end_dict[row['Gene']] = [row['Start_End']]
                #import pdb; pdb.set_trace()
            else:
                if row['Start_End'] not in gene_start_end_dict[row['Gene']]:
                    import pdb; pdb.set_trace()
                    target_genes_dict[row['Gene']][row['UTR']] ={}
                    target_genes_dict[row['Gene']][row['UTR']][row['Start_End']] ={}
                    target_genes_dict[row['Gene']][row['UTR']][row['Start_End']]['microRNA'] = row['microRNA']
                    target_genes_dict[row['Gene']][row['UTR']][row['Start_End']]['Seed'] = row['Seed']
                    target_genes_dict[row['Gene']][row['UTR']][row['Start_End']]['ddG'] = row['ddG']
                    gene_start_end_dict[row['Gene']].append(row['Start_End'])
                else:
                    import pdb; pdb.set_trace()
                    for utr , values in target_genes_dict[row['Gene']].items():
                        for key, value in values.items():
                            import pdb; pdb.set_trace()
                            if row['ddG'] < value['ddG'] :
                                pass
                #if [row['UTR'] in list(target_genes_dict[row['Gene']].keys()):

        '''
        target_csv_file = csv.reader(file, delimiter ='\t')
            for line in fh:
                line = line.strip()
                line_split = line.split('\t')
                try:
                    ddg_value = float(line_split[-1])
                except:
                    continue
                if min_energy >= ddg_value :
                    mi_rna = line_split[1]
                    if mi_rna not in target_genes_dict:
                        target_genes_dict[mi_rna] = {}
                    gene_notation = line_split[0]
                    if gene_notation not in target_genes_dict[mi_rna] :
                        target_genes_dict[mi_rna][gene_notation] ={}
                    target_genes_dict[mi_rna][gene_notation]['start'] = int(line_split[2])
                    target_genes_dict[mi_rna][gene_notation]['end'] = int(line_split[3])
                    target_genes_dict[mi_rna][gene_notation]['seed'] = line_split[4]
                    target_genes_dict[mi_rna][gene_notation]['ddG'] = ddg_value
        '''
    return target_genes_dict


if __name__ == '__main__' :
    if len(sys.argv) == 1:
        print('Usage summary_pita_results.py[arguments]')
        print('Type summary_pita_results.py --help for more information')
        exit(2)
    arguments = check_arg(sys.argv[1:])
    start_time = datetime.now()
    if not file_exists(arguments.reference):
        print ('3UTR reference genome does not exists\n')
        exit(2)
    if not directory_exists(arguments.dir):
        print ('directory with output results form PITA does not exists')
        exit(2)

    fetch_target_genes_data(arguments.dir, '_results.tab', -10.0)
