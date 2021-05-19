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
        num_match = 0
        not_match_num = 0
        match_file = os.path.join(self.__foldername, 'found_matches.tsv')
        no_match_file = os.path.join(self.__foldername, 'no_matches_found.txt')
        try:
            with open(match_file,'w') as fh_match, open(no_match_file, 'w') as fh_no_match :
                fh_match.write(self.__match_heading + '\n')
                fh_no_match.write('Organism\n')
                for elem in self.__buffer:
                    for key, values in elem.items():
                        if values != '':

                            # print('species : ' , key)
                            for value in values:
                                # print(value)
                                fh_match.write(str(key + '\t'+ value + '\n'))
                            num_match +=1
                        else:
                            fh_no_match.write(str(key + '\n'))
                            not_match_num +=1
        except:
            print('unable to create the output files')
            exit(2)
        print('\nSummary\n')
        print('Organisms found ', str(num_match))
        print('Number of Organisms with not match : ',str(not_match_num))



    def close(self):
        return super().close()


    def write(self, elem):  # resultado de la función BlastFunction
        self.__buffer.append(elem)
