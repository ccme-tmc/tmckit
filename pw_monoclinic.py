#!/usr/bin/env python
import list_utils as lu
import copy,sys
from py_utils import EmptyFormatter
from w2k_caseutil import *
from common_caseutil import KPointsT

def f_CreateNeighborVector(coords):
    '''
    Get all atom-atom vectors to determine whether two sets of atoms positions are in the same translation symmetry
    '''
    result = []
    for i in xrange(0,len(coords)):
        for j in xrange(i+1,len(coords)):
            vDiff = lu.f_List_Op_List(coords[i],"-",coords[j])
            result.append(f_MoveInPUC(vDiff))
            result.append(f_MoveInPUC(lu.f_List_Op_Scalar(vDiff,"*",-1)))
    result.sort()
    return result

def f_MoveInPUC(coord):
    '''
    Move a coordinate between [0,1] by add or sub 1
    '''
    coord2 = copy.copy(coord)
    for i in range(0,3):
        if ( coord[i] < 0):
            coord2[i] += 1
        elif ( coord[i] > 1):
            coord2[i] -= 1
    return coord2

def f_Matrix_compare(m1,m2,dErr=0.0001):
    '''
    Compare whether two matrixes are same
    '''
    for i in range(0,len(m1)):
        for j in range(0,len(m1[0])):
            if ( abs(m1[i][j]-m2[i][j]) > dErr):
                return False
    return True

def f_CompareTranslation(ref,new):
    '''
    Compare two sets of points and determine their relationship of translation, if exists
    :param ref: reference atoms, in internal coordinate
    '''
    ref.sort()
    for i in xrange(0,len(ref)):
#Try all translation vector constructed by ref[0] <- any new point
        vTrans = lu.f_List_Op_List(ref[0],"-",new[i])
        new2 = []
#Apply
        for j in xrange(0,len(ref)):
            new2.append(f_MoveInPUC(lu.f_List_Op_List(new[j],"+",vTrans)))

        new2.sort()
        
        print("Guess: %i to 0"% i)
        lu.f_Matrix_print(ref)
        lu.f_Matrix_print(new2)

        if ( f_Matrix_compare(ref,new2)):
            return new2

    #raise ValueError,"Cannot find suitable translation relationship between two set of points"

    return None

def w2k_GuessMonoclinicRotation(st1,st2):
    '''
    Try to guess the relationship between two monoclinc structure from cell vectors, which is represented by an unitary transformation AND a linear transformation ( O V2 = UT V1 U)
    This function is used to convert monoclinic structure before- and after- wien2k sgroup, to convert properties like k-points of before- to after-
    Note, as atom position is not involved in the guess, external procedure is required to verify the if the rotation follows the rule
    :param st1: cell vectors. each vector in one row.
    :param st2: the second vectors
    :return: rotation matrix U and linear transformation matrix O
    '''
    dErr = 0.0001
    #Correct nearly zero value to zero
    st1c = lu.f_Matrix_transpose(st1)
    st2c = lu.f_Matrix_transpose(st2)
    
    #lu.f_Matrix_print(st1c)
    #lu.f_Matrix_print(st2c)

    for i in xrange(0,3):
        for j in xrange(0,3):
            if ( abs(st1c[i][j]) < dErr):
                st1c[i][j] = 0.0
            if ( abs(st2c[i][j]) < dErr):
                st2c[i][j] = 0.0
    #Find two same axis
    listMap = []
    for i in xrange(0,3):
    #Only deal with two vectors exactly along xyz axis first
#This assume all common method to grab cell vectors from monoclinic cell will create two same vectors along axes
        bGoodVec = True
        for j in xrange(0,3): #skip that does not exactly mapped
            if ( i != j and st1c[j][i] != 0):
                 bGoodVec = False
                 break
        if ( not bGoodVec):
            continue
        for j in xrange(0,3):
            if ( abs(st1c[i][i]-st2c[j][j]) < dErr):
                listMap.append([i,j])
                break
    #If less than two is found then treat as an error
    if ( len(listMap) != 2):
        raise ValueError,"These two matrixes seems not two monoclinic matrix with some kind of relationship"


    #Create rotate matrix
    #print(listMap)
    listLeft1 = [0,1,2]
    listLeft2 = [0,1,2]
    for map in listMap:
        #R1[map[0]][map[1]] = 1
        listLeft1.remove(map[0]) 
        listLeft2.remove(map[1]) 
#Create the third map left
    listMap.append([listLeft1[0],listLeft2[0]])

#For every rotation mapping, we do not know whether it is positive or negtive, so guess for it
    
    Rbase = [[0,0,0],[0,0,0],[0,0,0]]

    listSign =[ [1,1,1],[1,1,-1],[1,-1,1],[-1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1],[-1,-1,-1]]
    #listSign = [[1,1,-1]]

    AInv = lu.f_List_inv3(st1c)
    #print(AInv)

    listResult = []
    for Sign in listSign:
    #set the remaining part to 1 or -1
        R1 = copy.deepcopy(Rbase)
        for aMap in zip(listMap,Sign):
            #print(aMap)
            R1[aMap[0][0]][aMap[0][1]] = aMap[1]
        #print("Rotate Matrix")
        #lu.f_Matrix_print(R1)
#Test which one is right
        O1 = lu.f_Matrix_dot(AInv,lu.f_Matrix_dot(lu.f_List_inv3(R1),lu.f_Matrix_dot(st2c,R1)))
#Detect whether O matrix include positive integer only
        #print("Transform Matrix")
        #lu.f_Matrix_print(O1)
        dSum1 = sum([sum([abs(x-abs(round(x,1))) for x in y]) for y in O1])

        if ( dSum1 < dErr):
            listResult.append([R1,O1])

    if ( len(listResult) == 0 ):
        raise ValueError,"Cannot find appropriate relationship between two monoclinic cell vectors"
    
#Return all the solutions
    return listResult


def w2k_ConvertMonoclinicKPoints(stFile1,stFile2,stKPtFile1,stKPtFile2):
    '''
    Convert k-points between two monoclinic structure
    '''
    w1 = Wien2K_Structure(stFile1)
    c1 = w1.GetCellRaw()
    mA = lu.f_Matrix_transpose(c1.PrimitiveCellVector)
    w2 = Wien2K_Structure(stFile2)
    c2 = w2.GetCellRaw()
    mB = lu.f_Matrix_transpose(c2.PrimitiveCellVector)
    listM = w2k_GuessMonoclinicRotation(c1.PrimitiveCellVector,c2.PrimitiveCellVector) 
    print("%i possible rotation found" % len(listM) )

#Test coord of atoms 
    nAtoms = len(w2.listAtom[0].InnerCoord) #The atoms should be considered in one translation verification ( there are no translation relationship between them, but what we only know is the non-equiv atom groups are same  in w1 and w2, but not single atoms)

    l1 = c1.GetAtomList(latt="prim",unit="bohr")[:nAtoms]
    l2 = c2.GetAtomList(latt="prim",unit="bohr")[:nAtoms]


    listC1 =[ x[1:4] for x in c1.listAtom[:nAtoms]] #Internal of c1 cell
    listDiff1 = f_CreateNeighborVector(listC1) 

    listT = []  #Store possible translation

    for M in listM: #enumerate all possible rotations

        print("Inner coord of B convert to A")
        mAinv = lu.f_List_inv3(mA)
        l2_A = [ lu.f_Matrix_transpose(lu.f_Matrix_dot(mAinv,lu.f_Matrix_dot(M[0],lu.f_Matrix_transpose([x[1:4]]))))[0]  for x in l2]

        lu.f_Matrix_print(listC1)
        mT =f_CompareTranslation(listC1,l2_A) 
        lu.f_Matrix_print(mT)
        listT.append(mT)

        if ( mT != None):
#Find an apporiate one
            KPt = w2k_GetKPointList(stFile1,stKPtFile1)
            KPt.ConvertUnit("cart",KPt.latt.PrimitiveCellVector)
            KPt2=copy.deepcopy(KPt)
            KPt2.listKPt = [ [x[0]] + lu.f_Matrix_transpose(lu.f_Matrix_dot(lu.f_List_inv3(M[0]),lu.f_Matrix_transpose([x[1:4]])))[0]+x[4:] for x in KPt.listKPt]
            KPt2.ConvertUnit("crystal",c2.PrimitiveCellVector)
            #KPt2.WriteToFile("test_sgroup.klist_band",True)
            w2k_WriteKPoints(KPt2,stKPtFile2)
            print("new k-points write into %s" % stKPtFile2)


def Main(ArgList):
    description = '''Convert k-points coordinate between two WIEN2k monoclinic structures. This is used to convert k-points if WIEN2k sgroup rotate monoclinic cell vectors from original one. In this case, the two structure files should be the one before sgroup that converted from cif, and after sgroup that used in calculation.
'''
    usage = "%prog -o KPointsFileInStructure2 Structure1File Structure2File KPointsFileInStructure1"
    
    parser = OptionParser(formatter=EmptyFormatter(),usage=usage,description=description)

    parser.add_option("-o",dest="KPointsFileName",default=None,help="The output k-points file name")

    (options,args) = parser.parse_args(ArgList)

    w2k_ConvertMonoclinicKPoints(args[1],args[2],args[3],options.KPointsFileName)

if __name__ == "__main__":
    Main(sys.argv);
