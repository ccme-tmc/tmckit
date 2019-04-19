#!/usr/bin/env python

# this Python module file contains functions related to Equation of States (EOS) fitting for a crystal 

import sys,os,shutil
import commands,string
from math import *
from numpy import *
from scipy.optimize import *

def f_EOS_Murnaghan(V,E0,V0,B0,B0P):
  return  E0 + (B0*V/B0P)*( (V0/V)**B0P/(B0P-1.0) + 1.0 ) - B0*V0/(B0P-1.0)

def f_EOS_BirchMurnaghan(V,E0,V0,B0,B0P):
  return  E0 + (9.0*V0*B0/16.0)*( ( (V0/V)**(2.0/3)-1.0 )**3*B0P + ( (V0/V)**(2.0/3) - 1.0 )**2*(6.0-4.0*(V0/V)**(2.0/3) ) )

def f_EOS_Harmonic(V,E0,V0,B0):
  return E0 + 0.5*B0*(V-V0)**2/V0

def f_EOS_Morse(R,E0,Re,De,a):
  return E0 + De*(1.0 - exp(-a*(R-Re)) )**2

def f_Resid_Murnaghan(p,etot,vol):
  V0 = p[0]
  E0 = p[1]
  B0 = p[2]
  B0P = p[3]
  resid = etot - f_EOS_Murnaghan(vol,E0,V0,B0,B0P)
  return resid

def f_Resid_BirchMurnaghan(p,etot,vol):
  V0 = p[0]
  E0 = p[1]
  B0 = p[2]
  B0P = p[3]
  resid = etot - f_EOS_BirchMurnaghan(vol,E0,V0,B0,B0P)
  return resid

def f_Resid_Morse(p,etot,R):
  (E0,Re,De,a) = (p[0],p[1],p[2],p[3])
  return etot - f_EOS_Morse(R,E0,Re,De,a) 

def f_EOS_Fitting(xval,etot,type,np_out=0,debug=False):

  # first do a polynomial fitting to get initial guess of E0, V0 and B0
  ndata = len(xval)

  xy_out = []
  if np_out > 0:
    x_min = min(xval) 
    x_max = max(xval)
    
    x0 = x_min - (x_max - x_min)*0.1
    x1 = x_max + (x_max - x_min)*0.1
    hx = (x1-x0)/(np_out-1)  
     
    for i in range(np_out):
      xy_out.append([x0+i*hx,0.0]) 
     

  coef = polyfit(asarray(xval),asarray(etot),2)
  c0 = coef[2]; c1 = coef[1]; c2 = coef[0]

  if type == 'M' or type == 'BM':  # EOS_Murnaghan
    E0 = c0-c1*c1/(4.0*c2)
    B0 = - c1
    V0 = - c1/(2.0*c2)
    B0P = 3.5
    print "Initial Guess from polyfit: \n C0=%12.6f, C1=%12.6f, C2=%12.6f" %(c0,c1,c2)
    print " E0=%12.6f, V0=%12.6f, B0=%12.6f" %(E0,V0,B0)

    p0 = array([ V0, E0, B0, B0P ])

    if type == 'M': 
      plsq = leastsq(f_Resid_Murnaghan, p0, args=(asarray(etot),asarray(xval)) )

      # calculate fitting error 
      err = 0.0
      for i in range(ndata):
        err += f_Resid_Murnaghan(plsq[0],etot[i],xval[i])**2
      err = sqrt( err/ndata )

    else:
      plsq = leastsq(f_Resid_BirchMurnaghan, p0, args=(asarray(etot),asarray(xval)))
      # calculate fitting error 
      err = 0.0
      for i in range(ndata):
        err += f_Resid_BirchMurnaghan(plsq[0],etot[i],xval[i])**2
      err = sqrt( err/ndata )

    par=plsq[0]
    V0 = par[0]; E0= par[1]; B0 = par[2]; B0P = par[3]
    print " Parameters from %s EOS fitting: E0=%16.6f, V0=%12.4f, B0=%12.6f and B0'=%12.6f" %(type,E0,V0,B0,B0P)
    if np_out > 0: 
      for i in range(np_out):
        V = xy_out[i][0]
        if type == 'M':
          xy_out[i][1]= f_EOS_Murnaghan(V,E0,V0,B0,B0P)
        else:
          xy_out[i][1]= f_EOS_BirchMurnaghan(V,E0,V0,B0,B0P)
    # 

  elif type == 'H':   # Harmonic EOS fitting 
    E0 = c0-c1*c1/(4.0*c2)
    x0 = -c1/(2.0*c2)
    k = 2*c2

    # calculate fitting error 
    err = 0.0 
    for i in range(ndata):
      err += ( etot[i] - (E0 + 0.5*k*(xval[i]-x0)**2) )**2
    err = sqrt( err/ndata )

    par = [ x0, E0,  k]


  elif type == 'Morse':

    Em = etot[0]
    Rm = xval[0]
    im = 0 
    for i in range(1,ndata):
      if etot[i] < Em:
        im = i
        Em = etot[i]
        Rm = xval[i] 

    if im == 0 or im == ndata-1:
      print "Warning: the minimum occurs at the boundary: Morse function is likely a poor fitting function!!!"
      if im == 0 : im += 1 
      if im == ndata-1: im -= 1
    
    # using three points that give lowest values to make a quadratic fitting 
    xtmp = array([xval[im-1],xval[im],xval[im+1]])
    ytmp = array([etot[im-1],etot[im],etot[im+1]])
    coef = polyfit(xtmp,ytmp,2)
    c0 = coef[2]; c1 = coef[1]; c2 = coef[0]

    E0 = c0-c1*c1/(4.0*c2)
    Re = - c1/(2.0*c2)
    k = 2*c2
    De = etot[ndata-1] - E0 
    a = sqrt(k/(2.0*De))

    if debug: 
      ofile=open("eosfit.dat",'w') 
      ofile.write("#Morse-type EOS:\n")
      ofile.write("# initial guess: E0= %12.6f, Re=%12.4f, De=%12.6f, a=%12.4f\n"%(E0,Re,De,a))

    p0 = array([ E0, Re, De, a]) 
    plsq = leastsq(f_Resid_Morse,p0,args=(asarray(etot),asarray(xval)))

    par = plsq[0]
    (E0, Re, De, a) = par 

    # calculate fitting error
    err = 0.0
    for i in range(ndata):
      err += f_Resid_Morse(par,etot[i],xval[i])**2
    err = sqrt( err/ndata )

    if debug:
      ofile.write("#  fit param.   : E0= %12.6f, Re=%12.4f, De=%12.6f, a=%12.4f\n"%(E0,Re,De,a))
      ofile.write("#  standard deviation=%12.6f\n"%(err))
      for i in range(ndata):
        x=xval[i]
        y=etot[i]
        yfit = f_EOS_Morse(x,E0,Re,De,a) 
        ofile.write("%16.6f %16.6f %16.6f\n"%(x,y,yfit))
      ofile.close() 

  err=err/(max(etot)-min(etot))
    
  return par,err,xy_out

