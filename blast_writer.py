from alemanpyutils.runner.writers.data_writer import DataWriter
import os

class BlastWriter(DataWriter):

    def __init__(self, foldername, match_heading):
        self.__foldername = foldername
        self.__buffer = []
        self.__match_heading = match_heading


    def open(self):
        return super().open()

    def pre(self):
        return super().pre()

    def post(self):
        num_match = [0,0]
        not_match_num = 0
        kind_matches =['direct','compl']
        match_file_direct = os.path.join(self.__foldername, 'found_matches_forward_direction.tsv')
        match_file_compl = os.path.join(self.__foldername, 'found_matches_complementary.tsv')
        no_match_file = os.path.join(self.__foldername, 'no_matches_found.txt')
        try:
            with open(match_file_direct,'w') as fh_match_direct, open(match_file_compl,'w') as fh_match_compl, open(no_match_file, 'w') as fh_no_match  :
                fh_match_direct.write(self.__match_heading + '\n')
                fh_match_compl.write(self.__match_heading + '\n')
                fh_no_match.write('Organism\n')
                for elem in self.__buffer:
                    for key in elem.keys():
                        found_any_match = 0
                        for k_match in kind_matches:
                            if elem[key][k_match] != '':
                                found_any_match = 1
                                for value in elem[key][k_match]:
                                    if k_match == 'direct':
                                    # print(value)
                                        fh_match_direct.write(str(key + '\t'+ value + '\n'))
                                        num_match[0] +=1
                                    else:
                                        fh_match_compl.write(str(key + '\t'+ value + '\n'))
                                        num_match[1] +=1
                        if found_any_match == 0:
                            fh_no_match.write(str(key + '\n'))
                            not_match_num +=1
        except:
            print('unable to create the output files')
            exit(2)
        print('\nSummary\n')
        print("Organisms found in 5' -- 3' direction: ", str(num_match[0]))
        print("Organisms found with complementary miRNA: ", str(num_match[1]))
        print('Number of Organisms which not match in any of the above: ',str(not_match_num))
        



    def close(self):
        return super().close()


    def write(self, elem):  # resultado de la función BlastFunction
        self.__buffer.append(elem)
