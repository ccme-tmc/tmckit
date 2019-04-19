#!/usr/bin/env python

import os
import commands


PBS_WriteDir = None #the directory to write file and run this script
PBS_RunDir = None #the working directory written into qsub file ( automatically stDir )  
PBS_Queue = "high" 

def f_CreatePBSScript(stCaseName,list_stCommand):
    '''
    Create a script that can be submit to pbs system , under specific dir
    :param stWorkDir: 
    :param stRunDir:  
    '''
    #stDir must be an absolute path
    stDir = PBS_WriteDir
    stRunDir = PBS_RunDir
    stQueue = PBS_Queue
    
    if ( stDir == None):
        stDir = os.getcwd()
    elif ( stDir[0] != "/" and stDir.find(":") == -1): # relative path
        stDir = os.path.join(os.getcwd(),stDir)
    
    if ( stRunDir == None):
        stRunDir = stDir
        
    fOut = open(os.path.join(stDir,"qsub-%s.sh" % stCaseName),"w")
    
    nProcess = 1
    if ( stQueue == "high"):
        nProcess = 24
    elif ( stQueue == "low"):
        nProcess = 16
        
    
    fOut.write('''#!/bin/sh -f
#PBS -o q-%s.out
#PBS -e q-%s.err
#PBS -l nodes=1:ppn=%d,walltime=100:00:00
#PBS -N %s
#PBS -q %s
#PBS -m ae

echo 'Job is running on node(s): '
cat $PBS_NODEFILE

export OMP_NUM_THREADS=1
export SCRATCH=%s

cd "%s"
'''
     % (stCaseName,stCaseName,nProcess,stCaseName,stQueue,os.environ["SCRATCH"],stRunDir)
    )
    for stCommand in list_stCommand:
        fOut.write("mpirun -hostfile $PBS_NODEFILE %s\n" % stCommand)
    
    fOut.close()
    
     


def f_GetEqualLineValue(stLine):
    '''
    Split a string like  "a = b" to b
    '''
    stPos = stLine.index("=")
    return stLine[:stPos].strip(),stLine[stPos+1:].strip()

class Node:
    '''
    Represent status of one node
    '''
    def __init__(self):
        self.stName = ""
        self.stState = "free"
        self.nProcess = 0
        self.stProperties = ""
        self.stType = "cluster"
        self.arJobs = {}
        self.stStatus = ""
    
    def ReadFromStringList(self,arLine,nIndex = 0):
        i = nIndex
        while(arLine[i] == ""):
            i += 1
        self.stName = arLine[i]
        
        i+= 1
        self.stState = f_GetEqualLineValue(arLine[i])[1]
        
        i += 1
        temp = f_GetEqualLineValue(arLine[i])[1]
        self.nProcess = int(temp)
        
        i += 1
        self.stProperties = f_GetEqualLineValue(arLine[i])[1]
        
        i += 1
        self.stType = f_GetEqualLineValue(arLine[i])[1]
        #Next : status, job
        #both can be skipped
        #Job can be empty, detect
        i += 1
        while ( True ):
            if ( arLine[i].strip() != ""):
                stName,stJobs =  f_GetEqualLineValue(arLine[i])
                if ( stName == "jobs"):
                    arJobs = [ [int(x.split("/")[0]),x.split("/")[1]] for x in stJobs.split(",") ]
                    self.arJobs = {}
                    for aJob in arJobs:
                        if ( not self.arJobs.has_key(aJob[1])):
                            self.arJobs[aJob[1]] = 1
                        else:
                            self.arJobs[aJob[1]] += 1
                    #stName,self.stStatus = f_GetEqualLineValue(arLine[i])
                else:
                    self.stStatus = stJobs
                i += 1
            else:
                break
        return i + 1
        
        
        

def f_ShowNodesDetail():
    status,output = commands.getstatusoutput("pbsnodes")

    if ( output[:20].find("pbsnodes") != -1):
        output = '''
pu4
     state = job-exclusive
     np = 16
     properties = seq
     ntype = cluster
     jobs = 0/563.plutonium, 1/563.plutonium, 2/563.plutonium, 3/563.plutonium, 4/563.plutonium, 5/563.plutonium, 6/563.plutonium, 7/563.plutonium, 8/563.plutonium, 9/563.plutonium, 10/563.plutonium, 11/563.plutonium, 12/563.plutonium, 13/563.plutonium, 14/563.plutonium, 15/563.plutonium
     status = opsys=linux,uname=Linux plutonium4 2.6.18-164.el5 #1 SMP Tue Aug 18 15:51:48 EDT 2009 x86_64,sessions=3078 5836 6308,nsessions=3,nusers=3,idletime=874910,totmem=131541516kb,availmem=129050876kb,physmem=66004360kb,ncpus=16,loadave=16.40,netload=5681372428585,state=free,jobs=563.plutonium,varattr=,rectime=1308819308

pu3
     state = job-exclusive
     np = 16
     properties = seq
     ntype = cluster
     jobs = 0/571.plutonium, 1/571.plutonium, 2/571.plutonium, 3/571.plutonium, 4/571.plutonium, 5/571.plutonium, 6/571.plutonium, 7/571.plutonium, 8/571.plutonium, 9/571.plutonium, 10/571.plutonium, 11/571.plutonium, 12/571.plutonium, 13/571.plutonium, 14/571.plutonium, 15/571.plutonium
     status = opsys=linux,uname=Linux plutonium3 2.6.18-164.el5 #1 SMP Tue Aug 18 15:51:48 EDT 2009 x86_64,sessions=6303 10237 10896 10850 23518,nsessions=5,nusers=4,idletime=67027,totmem=131541516kb,availmem=126797844kb,physmem=66004360kb,ncpus=16,loadave=16.64,netload=9821044465210,state=free,jobs=571.plutonium,varattr=,rectime=1308819313

node1
     state = free
     np = 24
     properties = para
     ntype = cluster
     jobs = 0/542.plutonium, 1/542.plutonium, 2/542.plutonium, 3/542.plutonium, 4/542.plutonium, 5/542.plutonium, 6/542.plutonium, 7/542.plutonium, 8/542.plutonium, 9/542.plutonium, 10/542.plutonium, 11/542.plutonium, 12/542.plutonium, 13/542.plutonium, 14/542.plutonium, 15/542.plutonium
     status = opsys=linux,uname=Linux node1 2.6.18-164.el5 #1 SMP Tue Aug 18 15:51:48 EDT 2009 x86_64,sessions=6440 18581,nsessions=2,nusers=2,idletime=1142660,totmem=61753140kb,availmem=61013468kb,physmem=49463424kb,ncpus=24,loadave=16.00,netload=4236896396747,state=free,jobs=542.plutonium,varattr=,rectime=1308819319

'''    
    
    arLine = [x.strip() for x in output.split("\n")]
    arNode = []
    i = 0
    while ( True):
        aNode = Node()
        i = aNode.ReadFromStringList(arLine, i)
        arNode.append(aNode)
        if ( i > len(arLine) - 4):
            break
    
    print(  "Node  Free  Status         Jobs")

    nLeft = 27 #The length of Node + Free + Status, the former two are the two most important parameter
    stLeft = "                           " 
    
    for aNode in arNode:
        stOut = "%-8s" % aNode.stName         
        stJobs = aNode.stState + "  "

        nCount = 0
        for stJobName,nJobUse in aNode.arJobs.iteritems():
            stJobs += stJobName+"*"+str(nJobUse)+","
            nCount += nJobUse
        if ( stJobs[-1] == ","):
            stJobs = stJobs[:-1]
        stOut += "%2d" %  (aNode.nProcess -nCount) + "  "
        stOut += stJobs
        # Wrap the line to avoid it appear in next line, use 80-width
        if ( len(stOut) > 79):
            listCutPos = [0] + [nLeft+x*(79-nLeft) for x in range(1,(len(stOut)-nLeft)/(79-nLeft)+1)] + [-1]
            stCut = []
            for i in range(0,len(listCutPos)-1):
                stCut.append(stOut[listCutPos[i]:listCutPos[i+1]])
            stOut = ( "\n"+stLeft).join(stCut)


        print(stOut)

if __name__ == "__main__":
    #f_ShowNodesDetail()
    try:
        f_ShowNodesDetail()
    except:
        print("Unable to query pbs status from server. Please use command 'pbsnodes'.")
    
     
        
