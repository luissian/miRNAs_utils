from alemanpyutils.function.function import Function
from Bio.Blast.Applications import NcbiblastnCommandline
import os
# from biopython import blast


class BlastFunction(Function):

    def __init__(self, params) -> None:
        self.__params = params
        parameter_list = params.replace(' ', '').split(',')[1:]
        self.__parameters_output = '\t'.join(parameter_list)

    def apply(self, elem):  #  elem = tupla que nos devuelve el reader

        query_fasta_file = elem[0]
        blast_db_name = elem[1]
        self.__params
        specie_name = os.path.basename(blast_db_name)
        print('procesing specie : ', specie_name)
        result = {specie_name:''}
        # resultado = blast(mirna_id, dir_id, params)
        #cline = NcbiblastnCommandline(db="/media/bioinfo/Massive_data/tesis_micro_RNA/parelizacion/test_genomas_ref_gastrointestinal/blastdb/Dysgonomonas_mossii_GCF_000376405/Dysgonomonas_mossii_GCF_000376405",task = 'blastn-short', evalue=0.001,  outfmt =6, num_threads=1, query=query_fasta_file)

        cline = NcbiblastnCommandline(db=blast_db_name, task='blastn-short', strand='minus', evalue=0.001,  outfmt =  self.__params ,num_threads=1, query=query_fasta_file)

        out, err = cline()
        out_lines = out.splitlines( )
        if len (out_lines) > 0 :
            matches_found = []
            for out_line in out_lines:
                line_split = out_line.split('\t')
                try:
                    if int(line_split[10]) > int(line_split[11]):
                        matches_found.append(out_line)
                except:
                    pass
            result = {specie_name:matches_found}
            #resultado = {specie_name:['hemos encontrado',str(len(out_lines)), 'matches']}


        return result
