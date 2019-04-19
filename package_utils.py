#!/usr/bin/env python


import os,sys,subprocess
from optparse import OptionParser

def f_Package_Detect(stDir=None,listFileName=None):
    '''
    Automatically determine which package may be used in specific folder or by list of filename .
    Note if two calculation cases are in the same folder, it may give wrong result.
    :param stDir: the directory should be detect. Default use current folder.
    :param listFileName : list of some file this package may used as input or output. Default None. If it is none or list with length=0, then directory information will be used, otherwise only this list is used
    :return: w2k,vasp,siesta,qe
    '''
    if ( listFileName== None ):
        listFileName = []
    if ( stDir == None and len(listFileName) == 0):
        stDir = os.getcwd()
        
    result = None

    f = open(os.devnull,'w')

    #Filename test
    if ( stDir != None):
        for root,dirs,files in os.walk(stDir):
            listFile = files
            break
    else:
        listFile = listFileName


    for stFileName in listFile:
        if ("w90" in stFileName):
            result = "w90"
            break
        if (stFileName.startswith("_gw") or stFileName.startswith("_GW")):
            result = "gap"
            break
        if ( "_DS" in stFileName):
            result = "abinit"
        if ( stFileName.upper() =="POSCAR"):
            result = "vasp"
            break
        else:
            stExt = os.path.splitext(stFileName.upper())[1]
            if ( stExt == ""):
                continue
            if ( stExt[0] == "."):
                stExt = stExt[1:]
#.struct is not used here as it is too common
            if ( stExt == "IN1" or stExt == "IN1C"):
                result = "w2k"
                break
            elif ( stExt == "FDF"):
                result = "siesta"
                break

    #File content test
    if ( result == None):
        for stFileName in listFile:
#qe input
            nExit = subprocess.call(["grep","-i","ATOMIC_SPECIES",stFileName],stdout=f,stderr=f)
            if ( nExit == 0):
                result = "qe"
                break
        #qe pw.x output
            nExit = subprocess.call(["grep","-i","PWSCF",stFileName],stdout=f,stderr=f) 
            if ( nExit == 0):
                result = "qe"
                break
#w2k struct
            nExit  =subprocess.call(["grep","-i","NONEQUIV",stFileName],stdout=f,stderr=f)
            if ( nExit == 0):
                result = "w2k"
                break
#yaeh band             
            nExit  =subprocess.call(["grep","-i","Special points were done with",stFileName],stdout=f,stderr=f)
            if ( nExit == 0):
                result = "yaeh"
                break
    
    f.close()

    return result

def uni_Package_Detect(stPackage="auto",listFileName=None):
    '''
    Wrap f_Package_Detect, automatically select package name if it is not determined
    :param stPackage: package name. if it is not "auto", then directly return
    :param listFileName : list of file that may be used by this package
    '''
    if ( stPackage == "auto"):
        stPackage = f_Package_Detect(listFileName=listFileName)
        if ( stPackage == None):
            raise ValueError,"Must specify package name"
        else:
            print("Detected package: %s" % stPackage)
    return stPackage

