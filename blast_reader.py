from alemanpyutils.runner.readers.data_reader import DataReader
import os



class BlastReader(DataReader):


    def __init__(self, all_mirnas_file:[], all_directories: []) -> None:
        self.__all_mirnas_file = all_mirnas_file
        self.__blastdb_main_dir = all_directories[0]
        self.__all_directories = all_directories[1]
        self.__count = 0
        self.__dir_index = 0


    def open(self):
        return super().open()

    def pre(self):
        return super().pre()

    def has_next(self):
        return self.__count < len(self.__all_directories)

    def post(self):
        return super().post()

    def close(self):
        return super().close()

    def next(self):
        next_species = os.path.join(self.__blastdb_main_dir,self.__all_directories[self.__dir_index],self.__all_directories[self.__dir_index] )
        next_elem = (
            self.__all_mirnas_file,
            next_species
            #self.__all_directories[self.__dir_index],
        )
        self.__dir_index += 1
        #self.__dir_index += 1
        self.__count += 1

        return next_elem
