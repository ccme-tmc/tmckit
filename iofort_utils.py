#!/usr/bin/env python
#This module helps to read input/output files of Fortran programs
import numpy as np
from struct import unpack

class FortranUnformattedFile(object):
    '''
    Deal with a Fortran unformatted file
    Note this file must be in ACCESS="SEQUENTIAL" mode, which is default
    some programs like vasp use "DIRECT" instead and not suitable
    '''
    def __init__(self,filename,mode='r'):
        '''
        :param filename: The filename to read
        '''
        self.filename = filename
        self.list_record = None
        if (mode == "r"):
            self.read()
    
    def read(self):
        self.f = open(self.filename,'rb')
        self.read_record()

    def read_record(self,b_print=False):
        '''
        Get the summary of records
        Also reset the pointer to records
        '''
        if (self.list_record is None):
            list_record = []
            f = self.f
            while True:
                pos = f.tell()
                t = f.read(4)
                if (t == ""):
                    break
                (l,) = unpack("i",t)
                f.seek(l,1)
                (l2,) = unpack("i",f.read(4))
                if (l != l2):
                    print("Found a bad record in the Fortran file %i %i" % (l,l2))
                list_record.append((pos,l))
            self.list_record = list_record
            
        if (b_print):
            for pos,l in self.list_record:
                print("%10i %10i" % (pos,l))

        self.i_record = -1

    def seek(self,n):
        '''
        Seek to the n-th record
        '''
        self.f.seek(self.list_record[n][0]+4,0)
        self.i_record = n

    def seeknext(self):
        '''
        Seek to the next record
        '''
        self.i_record += 1
        self.f.seek(self.list_record[self.i_record][0]+4,0)

    def readarray(self,dtype,shape,order='F'):
        '''
        Read an array from the file
        Note it is better to use numpy type with known bytes
        The order of the array is by default "F". The default of numpy is 'C'.
        This may change after arithmatic calculations...
        
        :param order: 'F' means read as-is, 'C' means read and convert to C order
        '''
        count = 1
        for n in shape:
            count *= n
        ar = np.reshape(np.fromfile(self.f,dtype,count),shape,order='F')
        if (order == 'C'):
            ar = np.ascontiguousarray(ar)
        return ar
    
    def close(self):
        '''
        '''
        self.f.close()
