#!/usr/bin/python3


import multiprocessing
import os  

from glob import glob
from tirAnalysis.irbisCSV2tiff import (irbisCSV2data,
                                       irbisData2Tiff)


def process_file(filename, root='', limits=(None, None), cmap='gray'):
    '''
    Wrapper function that declares the file-processing steps
    '''
    outpath = '../tiffs/'
    if root is not '':
        outpath += root + '/'
    if not os.path.isdir(outpath):  
        os.mkdir(outpath)  

    data, md = irbisCSV2data(filename, sep=',')  
    nu_fname = outpath + filename[:-4]  
    print(filename + '\t', end='')  
    irbisData2Tiff(nu_fname, data, cmap=cmap, limits=limits)

    return


def multiprocessor_irbis2tiff(ROOT='', function='process_file'):
    '''
    Run a csv-to-image producing function on multiple csv files, each on its 
    own processor.  
    '''
    p = multiprocessing.Pool() 
    for filename in sorted(glob(ROOT + '*csv')):  
        p.apply_async(function, [filename, ROOT, (15, 45)])  
    p.close()      
    p.join()

    return


if __name__ == '__main__':
    ROOT = ''
    multiprocessor_irbis2tiff(ROOT)
