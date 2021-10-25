#!/usr/bin/python3

""" readCSVfiles 
"""

import numpy as np


def irbisCSV2tiff(fname):
    f = open(fname, 'r', encoding='utf-8', errors='ignore')

    found = False
    data  = []
    nread = 0
    metadata = {}
    mdread = 0

    ##--------------------READ METADATA--------------------##
    key   = None
    while key != '':
        # Remove leading and trailing whitespace, as well as any newline
        # characters.  Split the remaining output.
        key, *params = map(str.strip, f.readline().split('='))
        if len(params) > 0: 
            params = params[0].split(';')
        if mdread > 0 and key != '':
            metadata[key] = params
        mdread += 1

    ##--------------------READ THE DATA--------------------## 
    for line in f:
        if found:
            line = line.strip().split('\t')
            data.append(line)
            nread += 1  # should equal 480 if the code completes successfully
        else:
            if line.strip() == '[Data]': found = True

    f.close()

    return np.float64(data), metadata


##------------------------------------MAIN------------------------------------##
if __name__ == "__main__":
    from tifffile import imsave

    for num in range(1, 1401):
        fname = './AA103151_%04d.csv' % num
        # cut the .csv extention and add .tif
        foutn = '{0}.tif'.format(fname[:-4])  
        data, _ = irbisCSV2tiff(fname)
        imsave(foutn, data)

