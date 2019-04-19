#!/usr/bin/env python
import sys,os,shutil
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
  return acos(1.0*f_List_dot3(v1,v2)/f_List_norm3(v1)/f_List_norm3(v2))

def f_List_norm3(v):
  '''
  the norm of a 3-d vector
  '''
  f_List_Check(v)
  return v[0]**2 + v[1]**2 + v[2]**2
    

def f_List_Check(v):
  '''
  Check if a parameter is 3-d vector, raise error if it is not
  '''
  err = 0
  if ( (not isinstance(v,list)) ):
    err = 1
  elif ( len(v) != 3):
    err = 1
  if ( err == 1):
    print("Unsupported value:",v)
    raise TypeError,"Only 3-dimensional vector is supported"


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

