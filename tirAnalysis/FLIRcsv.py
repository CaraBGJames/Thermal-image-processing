import numpy as np


def readFLIRE4csv(filename, sep=';'):
    '''
    Returns an (numpy) array containing the radiometric data from FLIR E4 csv 
    files.

    Parameters:
    -----------
    filename (str)
        Absolute or relative path of the file to read
    sep (str)
        Field separator.  Default is semicolon ';'


    Returns:
    --------
    data (np.array)
        Floating point (radiometric) image

    '''
    readData = False

    data  = []
    count = 0

    with open(filename, 'r') as f:
        for ind, line in enumerate(f):
            if line[:7] == 'Image 1':
                readData = True
            if readData:
                line = line.strip().split(sep)
                if line[-1] == '':
                    line = line[:-1]
                    readData = False # stop reading at this point
                if line == '':
                    break
                else:
                    data.append(line[1:])
                    count += 1

    data = np.float64(data)
    print("%d lines read" % count)
    return data
