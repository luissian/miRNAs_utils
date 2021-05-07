from alemanpyutils.runner.writers.data_writer import DataWriter


class BlastWriter(DataWriter):

    def __init__(self, filename):
        self.__filename = filename
        self.__buffer = []


    def open(self):
        return super().open()

    def pre(self):
        return super().pre()

    def post(self):
        count = 0
        for elem in self.__buffer:
            for key, values in elem.items():
                if values != '':
                    print('species : ' , key)
                    for value in values:
                        print(value)
                    count +=1
        for elem in self.__buffer:
            for key, values in elem.items():
                if values == '':
                    print('Not match for species : ', key)
        print('Organismos que encuentran en blast ', str(count))


    def close(self):
        return super().close()


    def write(self, elem):  # resultado de la función BlastFunction
        self.__buffer.append(elem)
