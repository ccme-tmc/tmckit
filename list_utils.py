#!/usr/bin/env python
import sys,os,shutil,copy
from math import *


# This defines some home-made functions to manipulate lists to partially 
# replace functions in the numpy. This should be used only when the dimension of 
# the lists involved is very small 

def f_List_inv3(m):
  """
  Invert a 3x3 matrix represented by a row-wise list 
  i.e. m[i][:] refers to i-th volume in the matrix  
  """

  (a11,a21,a31) = (m[0][0], m[0][1], m[0][2])
  (a12,a22,a32) = (m[1][0], m[1][1], m[1][2])
  (a13,a23,a33) = (m[2][0], m[2][1], m[2][2])

  denom = a11*a22*a33 + a21*a32*a13 + a31*a12*a23 - a31*a22*a13 - a21*a12*a33 - a11*a32*a23
  denom *= 1.0  #force float type
  
  minv = []
  minv.append( [ a22*a33 - a23*a32,  a23*a31 - a21*a33, a21*a32 - a22*a31 ] )

  minv.append( [ a13*a32 - a12*a33,  a11*a33 - a13*a31, a12*a31 - a11*a32 ] )

  minv.append( [ a12*a23 - a13*a22,  a13*a21 - a11*a23, a11*a22 - a12*a21 ] )

  for i in range(3):
    for j in range(3):
      minv[i][j] /= denom 

  return minv 

def f_List_mv3(m,v):
    """
    the product between a 3x3 matrix, in the row-wise list format, and a vector
    """
    mv=[0.0,0.0,0.0]
    for i in range(3):
        for j in range(3):
            mv[i] += m[j][i]*v[j]  
    return mv 

def f_List_dot3(v1,v2):
    '''
    the dot product of 2 vectors
    '''
    f_List_Check(v1)
    f_List_Check(v2)
    return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

def f_List_angle3(v1,v2):
    '''
    the angle between 2 vectors, unit in rad  
    '''  
    f_List_Check(v1)
    f_List_Check(v2)
    d_acos = 1.0*f_List_dot3(v1,v2)/f_List_norm3(v1)/f_List_norm3(v2)
    if ( d_acos > 1):
        d_acos = 1
    elif ( d_acos < -1):
        d_acos = -1
    return acos(d_acos)

def f_List_norm3(v):
    '''
    the norm of a 3-d vector
    '''
    f_List_Check(v)
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)

def f_List_cross3(v1,v2):
    '''
    the cross product of 2 3-d vector
    '''
    f_List_Check(v1)
    f_List_Check(v2)
    return [v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0]]
  
def f_List_norm(v):
    '''
    the norm of any dimensional vector
    '''
    result = 0.0
    for x in v:
        result  = result + x*x
    result = sqrt(result)
    return result

def f_Matrix_print(m1):
    '''
    Print a matrix in human-readable format
    '''
    if (m1 == None):
        print("Empty matrix")
        return
    if(len(m1) == 0):
        print("Empty matrix")
        return
    elif ( len(m1[0]) == 0):
        print("Empty matrix")
        return

    stFormat = "|"+" ".join(["%8.4f" for x in m1[0]])+" |"
    stHead = stFormat.replace("%8.4f","        ")
    print(stHead)
    for row in m1:
        print(stFormat % tuple([ float(x) for x in row]))
    print(stHead)

def f_Matrix_dot(m1,m2):
    '''
    Dot product of two arrays
    '''
    l1_1 = len(m1)
    l1_2 = len(m1[0])
    l2_1 = len(m2)
    l2_2 = len(m2[0])
    if ( l1_2 != l2_1):
        raise ValueError,"The column count of matrix 1 (%d) must equal to row count of matrix 2 ( %d )" % (l1_2,l2_1)
    m3 = f_Matrix_zeros(l1_1,l2_2)
    for i in range(0,l1_1):
        for j in range(0,l2_2):
            for k in range(0,l1_2):
                m3[i][j] += m1[i][k]*m2[k][j]
    return m3

def f_Matrix_dimension(m):
    '''
    Check if a variable is matrix and return dimension x,y
    '''
    if ( isinstance(m,list)):
        if ( isinstance(m[0],list)):
            return len(m),len(m[0])
    
    raise TypeError,"Input variable is not a list-list like matrix"


def f_Matrix_transpose(m):
    '''
    the transpose of matrix
    '''

    l1 = len(m) # row count
    if ( isinstance(m[0],list)):
        l2 = len(m[0]) # column count
        m2 = m
    else:
        raise TypeError,"Input variable is not a list-list like matrix"
    #else:
    #    l1 = 1
    #    l2 = len(m)
    #    m2 = [m]
    #m2 is new input matrix
    
    
    m3 = f_Matrix_zeros(l2,l1)
    for i in range(0,l2):
        for j in range(0,l1):
            m3[i][j] = m2[j][i]            
            
    #if ( l2 == 1): # single column -> single row,
    #    m3 = m3[0]
        
    return m3
            
def f_Matrix_zeros(l1,l2):
    '''
    Create a l1 * l2 array filled with 0
    '''
    return [[0.0 for x in range(0,l2)] for y in range(0,l1)]

def f_Matrix_Op_Scalar(m,op,s):
    l1 = len(m)
    l2 = len(m[0])
    m2 = copy.deepcopy(m)
    for i in range(0,l1):
        for j in range(len(m[0])):
            if op == '+': m2[i][j] += s
            elif op == '*': m2[i][j] *= s 
            elif op == '-': m2[i][j] -= s
            elif op == '/': m2[i][j] /= s
            elif op == '**': m2[i][j] = m[i][j]**s 
            else: 
              raise ValueError,"unsupported operator in f_Matrix_Op_Scalar"
    return m2      

def f_List_permute(ar,nCount):
    '''
    Premute an array. [0,1,2] will be [2,0,1] if permute 1 
    '''
    nMax = len(ar)
    if ( nMax <= 1 or nCount == 0 ):
        return ar

    ar2 = copy.copy(ar)
    
    for i in range(0,nMax):
        nIndex = (i - nCount) % nMax
        ar2[i] = ar[nIndex]
    
    return ar2

def f_List_normalize(v):
    '''
    normalize a vector make norm(vec) = 1
    '''
    f_List_Check(v)
    n = f_List_norm(v)*1.0
    result = []
    for value in v:
        result.append(value/n)
        
    return result

def f_List_IsParallel(v1,v2):
    '''
    Determine whether two vectors are parallel
    '''
    f_List_Check(v1)
    f_List_Check(v2)
    dErr = 0.0000001
    if ( abs(f_List_norm3(f_List_cross3(v1,v2))) < dErr  ):
        return True
    else:
        return False

def f_Matrix_det(m):
    '''
    Get determinant of a [n,n] matrix
    '''
    n1,n2 = f_Matrix_dimension(m)
    if ( n1 != n2 ):
        raise TypeError,"Determinat must be sqaure!"
    if ( n1 == 1 ):
        return m[0][0]
    else:
        result = 0
        for i in range(0,n1):
            m2 = copy.deepcopy(m[1:])
            for j in range(0,n1-1):
                m2[j][i:i+1] = []
            result += (m[0][i] * f_Matrix_det(m2)) * ( -1 ) ** i
    return result


def f_List_Check(v):
  '''
  Check if a parameter is 3-d vector, raise error if it is not
  '''
  err = 0
  if ( not ( isinstance(v,list) or isinstance(v,tuple)  )):
    err = 1
  elif ( len(v) != 3):
    err = 1
  if ( err == 1):
    print("Unsupported value:",v)
    raise TypeError,"Only 3-element vector is supported"


def f_Str_Included(str,word):
  try:
    id = str.split().index(word)
    return True
  except:
    return False

def f_Str_New(s0,n):
  """
  Create a string with repeated s0
  """
  s=''
  for i in range(n): s += s0
  return s

def f_List_New(a,n):
  """ 
  Create a list with idential elements a
  """
  l = []
  for i in range(n):
    l.append(a)
  return l

def f_List_Op_Scalar(l,op,s):
  lnew=l[:]
  for i in range(len(lnew)):
    if op == '+': lnew[i] += s
    elif op == '*': lnew[i] *= s 
    elif op == '-': lnew[i] -= s
    elif op == '/': lnew[i] /= s
    elif op == '**': lnew[i] = l[i]**s 
    else: 
      print "ERROR in f_List_Op_Scalar: unsupported option!"
      sys.exit(1) 
  return lnew    

def f_List_Op_List(l1,op,l2):
  lnew=l1[:]
  n = len(lnew) 
  if op == '+': 
    for i in range(n): lnew[i] += l2[i]
  elif op == '-': 
    for i in range(n): lnew[i] -= l2[i]
  elif op == '*': 
    for i in range(n): lnew[i] *= l2[i]
  elif op == '/': 
    for i in range(n): lnew[i] /= l2[i]
  elif op == '.': 
    sum =0
    for i in range(n): sum += l1[i]*l2[i]
    return sum 
  else: 
    print "ERROR in f_List_Op_List: unsupported option %s!"%(op)
    sys.exit(1)
  return lnew

def f_List_Included(elem_new,elems,tol=1.e-10):
  """
  Check whether elem_new is already included in the list "elems"
  """
  dist_min = 1.e20 
  i_min = -1
  for i in range(len(elems)): 
    dist = 0.0 
    for ix in range(len(elem_new)): dist += (elem_new[ix] - elems[i][ix])**2
    dist = sqrt(dist)
    if dist < dist_min: 
      dist_min = dist 
      i_min = i

  if dist_min  > tol: i_min = -1

  return i_min 

