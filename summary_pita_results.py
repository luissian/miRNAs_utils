import sys, os, glob, re
import pandas as pd
pd.options.mode.chained_assignment = None
import argparse
from Bio import SeqIO
from datetime import datetime
from utils.generic_functions import *

'''

usage example
python3 ~/tesis_code/miRNAs_utils/summary_pita_results.py  -d pita_results_15-jan -r /home/lchapado/bio_programs/PITA/3UTRhuman.fa -p aligment -e -10 -m 2

ath_miR165a_5p_pita_results.tab
requires to install xlsxwriter
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

    parser.add_argument('-m','--min', required = True, help = 'minimun number of target genes')

    return parser.parse_args()

def convert_pita_dict_to_panda (target_genes_dict):
    '''
    Description:
        The function get the dicitonnary and conver it to panda dataframe.
    Input:
        target_genes_dict   # directory to convert
    Return:
        target_genes_dict
    Gene    Annotation  Start   End     microRNA    Seed    ddG
    UTR	microRNA	Start	End	Seed	Loop	dGduplex	dG5	dG3	dG0	dG1	dGopen	ddG
    ENST00000310015|SP3	ath_miR165a-5p	6872	6864	8:0:1	0	-21.1	-7.1	-14	-29.49	-27.61	-1.87	-19.22
    ENST00000320285|AGPAT4	ath_miR165a-5p	4781	4773	8:0:1	0	-24.52	-9.3	-15.22	-41.37	-35.52	-5.85	-18.66
    '''
    heading = {'Gene':[], 'Annotation':[], 'Start':[], 'End':[], 'microRNA':[],  'Seed':[], 'ddG':[]}
    converted_df = pd.DataFrame(heading)
    for gene, utrs in target_genes_dict.items():
        for utr in utrs.keys():
            for index in range(len(utrs[utr]['data'])):
                start, end = utrs[utr]['data'][index]['Start_End'].split('_')
                new_row = {'Gene':gene, 'Annotation':utr, 'Start':start, 'End':end, 'microRNA':utrs[utr]['data'][index]['microRNA'],
                            'Seed':utrs[utr]['data'][index]['Seed'], 'ddG':utrs[utr]['data'][index]['ddG']}
                converted_df = converted_df.append(new_row, ignore_index = True)

    return converted_df



def fetch_target_genes_data(files_to_process,  min_energy):
    '''
    Description:
        The function get the target genes form the PITA output files and reduce the
        results whith the ones that have less energy value as indicated in min_energy
        and it stores in target_genes_dict
    Input:
        files_to_process   # output PITA directory
        min_energy      # minimun value to filter the target genes
    Return:
        target_genes_dict
    UTR	microRNA	Start	End	Seed	Loop	dGduplex	dG5	dG3	dG0	dG1	dGopen	ddG
    ENST00000310015|SP3	ath_miR165a-5p	6872	6864	8:0:1	0	-21.1	-7.1	-14	-29.49	-27.61	-1.87	-19.22
    ENST00000320285|AGPAT4	ath_miR165a-5p	4781	4773	8:0:1	0	-24.52	-9.3	-15.22	-41.37	-35.52	-5.85	-18.66
    '''
    target_genes_dict = {}
    gene_start_end_dict = {}

    for file in files_to_process:
        panda_df = pd.read_csv(file, sep ='\t', header=0)
        selected_df = panda_df[panda_df['ddG'] <= min_energy]

        selected_df[['UTR','Gene']] = selected_df['UTR'].str.split('|',expand = True)
        #selected_df[['UTR','Gene']] = selected_df.loc[:,['UTR']].str.split('|',expand = True)
        selected_df['Start_End'] = selected_df['Start'].apply(str).str.cat(selected_df['End'].apply(str),sep = '_')
        for index, row in selected_df.iterrows():
            if not row['Gene'] in target_genes_dict:
                target_genes_dict[row['Gene']] = {}
                gene_start_end_dict[row['Gene']] = []

            if row['Start_End'] not in gene_start_end_dict[row['Gene']]:
                if not  row['UTR'] in target_genes_dict[row['Gene']] :
                    target_genes_dict[row['Gene']][row['UTR']] ={}
                    target_genes_dict[row['Gene']][row['UTR']]['positions']= []
                    target_genes_dict[row['Gene']][row['UTR']]['data'] = []
                data = {}
                data['microRNA'] = row['microRNA']
                data['Start_End'] = row['Start_End']
                data['Seed'] = row['Seed']
                data['ddG'] = row['ddG']

                target_genes_dict[row['Gene']][row['UTR']]['data'].append(data)
                target_genes_dict[row['Gene']][row['UTR']]['positions'].append(row['Start_End'])
                gene_start_end_dict[row['Gene']].append(row['Start_End'])
            else:
                # find the existing value
                found = False
                for utr in target_genes_dict[row['Gene']].keys():
                    if not row['Start_End'] in target_genes_dict[row['Gene']][utr]['positions']:
                        continue

                    for index_data in range(len(target_genes_dict[row['Gene']][utr])):
                        if not row['Start_End'] == target_genes_dict[row['Gene']][utr]['data'][index_data]['Start_End'] and row['microRNA'] == target_genes_dict[row['Gene']][utr]['data'][index_data]['microRNA'] :
                            continue
                        found = True
                        if row['ddG'] < target_genes_dict[row['Gene']][utr]['data'][index_data]['ddG']:
                            target_genes_dict[row['Gene']][utr]['data'][index_data]['ddG'] = row['ddG']
                            target_genes_dict[row['Gene']][utr]['data'][index_data]['Seed'] = row['Seed']
                        break
                    if not found:
                        data = {}
                        data['microRNA'] = row['microRNA']
                        data['Start_End'] = row['Start_End']
                        data['Seed'] = row['Seed']
                        data['ddG'] = row['ddG']
                        target_genes_dict[row['Gene']][utr]['data'].append(data)
                        break

    return target_genes_dict

def filter_min_target_genes(target_genes_dict, min_number_target_gene):
    '''
    Description:
        The function filter the target genes, returning  a dataframe with only the ones that the matches
        is equal or bigger to min_number_target_gene
    Input:
        files_to_process   # output PITA directory
        min_energy      # minimun value to filter the target genes
    Return:
        converted_df
    '''
    converted_df ={}
    heading = {'Gene':[], 'Annotation':[],'microRNA':[], 'Number_of_sites':pd.Series([], dtype='int'), 'Start_End':[]}
    converted_df = pd.DataFrame(heading)
    for gene, utrs in target_genes_dict.items():

        number_found = 0
        if len(utrs) >= min_number_target_gene:
            number_found = len(utrs)
        else:
            for utr in utrs.keys():
                number_found += len(utrs[utr]['positions'])
        if number_found >= min_number_target_gene:
            for utr in utrs.keys():
                n_sites = len(utrs[utr]['positions'])
                new_row = {'Gene':gene, 'Annotation':utr, 'microRNA':utrs[utr]['data'][0]['microRNA'], 'Number_of_sites': n_sites, 'Start_End':','.join(utrs[utr]['positions'])}
                converted_df = converted_df.append(new_row, ignore_index = True)

    return converted_df

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
    try:
        score = float(arguments.score)
    except:
        print('score must be a number')
        exit(2)
    try:
        min_number_target_gene = int(arguments.min)
    except:
        print('minimun number of target genes must be an integer')
        exit(2)

    filter_name = '_results.tab'
    out_file_name = 'prueba_multi.xlsx'
    excel_file = pd.ExcelWriter(out_file_name, engine = 'xlsxwriter')
    filter_files = os.path.join(arguments.dir,'*' + filter_name)
    files_to_process =glob.glob(filter_files)
    for file in files_to_process:
        print('Processing ',file)
        target_genes_dict = fetch_target_genes_data([file], score)
        filter_min_target_df = filter_min_target_genes(target_genes_dict, min_number_target_gene)
        gene_target_df = convert_pita_dict_to_panda(target_genes_dict)
        sheet_name = os.path.basename(file).replace(filter_name,'')
        save_dataframe_to_excel_multiple_sheet(gene_target_df,excel_file,sheet_name)
        save_dataframe_to_excel_multiple_sheet(filter_min_target_df,excel_file,'min_target_' + sheet_name)
    # add all miRNAs files
    if len(files_to_process) >1:
        print('Processing all miRNAs')
        target_genes_dict = fetch_target_genes_data(files_to_process, score)
        gene_target_df = convert_pita_dict_to_panda(target_genes_dict)
        save_dataframe_to_excel_multiple_sheet(gene_target_df,excel_file,'all_miRNAs')
    excel_file.save()
