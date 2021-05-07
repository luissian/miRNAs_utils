
from alemanpyutils.runner.runner_config import RunnerConfig
from alemanpyutils.runner.readers.reader import Reader
#from alemanpyutils.runner.readers.random_string_data_reader import RandomStringDataReader
from alemanpyutils.runner.writers.writer import Writer
#from alemanpyutils.runner.writers.standard_output_writer import StandardOutputWriter
from alemanpyutils.runner.processors.processor import Processor
#from alemanpyutils.function.to_upper_case_function import ToUpperCaseFunction
from alemanpyutils.runner.runner import Runner
from blast_reader import BlastReader
from blast_writer import BlastWriter
from blast_function import BlastFunction
import sys, os, argparse, pdb
from Bio import SeqIO

'''
Este es mi fichero de prueba para probar la realización de hacer blast para
las especies de la microbiota.
creamos la clase BlastReader que es donde definiremos la lista de los miRNAS y
los directorios de Blastn

Creamos la clase BlastFunction que es donde esta el grueso de la función para hacer blast
y luego coger la información

Creamos la clase BlastWriter que es la encargada de crear el fichero juntando toda
la información de la clase BlastFuncion

Exampld:
python mytest.py  -b /media/bioinfo/Extreme\ SSD/Tesis_uam/genomas_ref_gastrointestinal/blastdb/ -m /media/bioinfo/Massive_data/Reference_genomes/miRNAs_mirbase/mature_hsa.fasta

'''

def check_arg (args=None) :
    '''
    The function is used for parsing the input parameters form the command line
    using the standard python package argparse. The package itself is handling
    the validation and the return errors messages
    Input:
        args    # Contains the arguments from the command line
    Return:
        parser.parse_args()     # The variable contains the valid parameters
    '''
    parser = argparse.ArgumentParser(prog = 'mytest.py',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description= 'mytest is a test to implement paralelization ')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')
    parser.add_argument('-b', '--blastdbDir' , required = True,
                        help='blastdb folder with all organism index')
    parser.add_argument('-m', '--mirnafile', required =True,
                        help='miRNA file')
    return parser.parse_args()

def check_if_file_exists(file_to_check):
    '''
    Description:
        The function will check if the file exists
    Input:
        file_to_check       # file name
    Return:
        True if exists , else False
    '''
    if os.path.isfile(file_to_check):
        return True
    else:
        return False

def check_if_folder_exists(folder_to_check):
    '''
    Description:
        The function will check if the file exists
    Input:
        folder_to_check       # folder name
    Return:
        True if exists , else False
    '''
    if os.path.isdir(folder_to_check):
        return True
    else:
        return False

def read_mirna_file(mi_rna_file):
    '''
    Description:
        The function reads the miRNA file and return a tupla with miRNA_name and miRNA_sequence
    Input:
        mi_rna_file       # miRNA file name
    Return:
        mi_rna_list
    '''
    mi_rna_list = []
    for seq_record in SeqIO.parse(mi_rna_file, 'fasta'):
        mi_rna_list.append([str(seq_record.id), str(seq_record.seq)])

    return mi_rna_list


if __name__ == '__main__':
    if len (sys.argv) == 1 :
        print('Usage: find_validated_genes.py [ARGUMENTS] ')
        print('Try  find_validated_genes.py --help for more information.')
        exit(2)
    arguments = check_arg(sys.argv[1:])
    if ' ' in arguments.blastdbDir:
        print('Blastdb folder cannot have spaces')
        exit(2)
    if not check_if_folder_exists(arguments.blastdbDir):
        print('Blastdb folder does not exists\n')
        exit(2)
    if not check_if_file_exists(arguments.mirnafile):
        print('miRNA file does not exists\n')
        exit(2)
    blastdb_folder = os.listdir(arguments.blastdbDir)

    ## necesitamos el fichero no la secuencia
    # mi_rnas = read_mirna_file(arguments.mirnafile)
    blast_parameters = '"6 , qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend , sseq , qseq"'
    rc = RunnerConfig(batch_size=4, num_processors = 4)
    reader = Reader(data_reader = BlastReader(arguments.mirnafile, (arguments.blastdbDir, blastdb_folder)),runner_config = rc)
    writer = Writer(data_writer=BlastWriter('test_file'), runner_config= rc)

    processors = []
    print('numero especies' , len(blastdb_folder))
    for i in range(0,rc.num_processors):
        processors.append(Processor(reader, BlastFunction(blast_parameters), writer,rc))

    runner = Runner(reader,writer,processors,rc)
    runner.launch()
