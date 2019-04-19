#!/usr/bin/env python
import sys,os,shutil
import commands,string
from math import *
from subprocess import *
from scipy.optimize import *
from eos_fitting import *
from numpy import *
from list_utils import *

Pi = 3.14159265359
GoldenRatio = (sqrt(5.0)-1)/2.0

opt_restart_file="opt.restart"

def wrap_function(function, args): 
  def function_wrapper(x):
    return function(x, *args)
  return function_wrapper

def f_Read_Restart():
  # enable restart mode
  ifile = open(opt_restart_file,'r')
  lines = ifile.readlines()
  ifile.close()

  line = lines[len(lines)-1]
  s=line.split()
  ndim = len(s)-1
  x=[]
  for i in range(ndim): x.append(float(s[i]))
  fx = float(s[ndim])
  return x,fx

# write out results of the current cycyle to be used for restarting 
def f_Write_Restart(x,fx):
  ndim = len(x)
  ofile = open(opt_restart_file,'w')
  for i in range(ndim): ofile.write(" %24.12e " %(x[i]))
  ofile.write(" %24.12e \n" %(fx))
  ofile.close()

def f_Min1D_GoldenSearch(func,a,b,args=(),nmax=100,xtol=1.e-6,ftol=1.e-6,debug=False,save_all=False):
  """
  Search the minimum of a 1-D function by the Golden bracketing method 
  """
  func = wrap_function(func,args)
  r1 = GoldenRatio
  r2 = r1*r1
  if a > b: t = a; a = b; b= t
     
  h = b-a 
  ya = func(a)
  yb = func(b)

  c = a + r2*h
  yc = func(c)
  # check whether the function is really bracketed 
  if yc > ya or yc > yb:
    print "ERROR in f_Min1D_GoldenSearch: the function is not bracked in [%12.6f,%12.6f]"%(a,b)
    sys.exit(1)

  d = a + r1*h 
  yd = func(d)

  if save_all:
    save_data=[(a,ya),(c,yc),(d,yd),(b,yb)]

  if yc < yd:
    x0 = c; y0 = yc
  else:
    x0 = d; y0 = yd

  i = 1 
  while i < nmax: 
    y0_old = y0 

    if yc < yd: 
      b = d; yb = yd
      d = c; yd = yc  
      h = b - a
      c = a + r2*h 
      yc = func(c)
      x0 = c; y0 = yc

    else:
      a = c; ya = yc
      c = d; yc = yd
      h = b-a
      d = a + r1*h
      yd = func(d)
      x0 = d; y0 = yd

    if save_all:
      save_data.append((x0,y0))

    i += 1
    if debug:
      print "i=%5d, a=%12.6f, f(a)=%16.6f; b=%12.6f, f(b)=%16.6f" %(i,a,ya,b,yb)

    if (b-a) < xtol or abs(y0-y0_old) < ftol: break  

  if i == nmax:
    print "WARNING: convergence not reached after %d iterations" %(nmax)  

  if save_all:
    return x0,y0,save_data
  else:
    return x0,y0

def f_Min1D_Bracket(func,ax,bx,args=(),xmin=-1.1e10,xmax=1.1e10,glimit=10.0,debug=False,save_all=False):
  """
  This function is a modfied version of mnbrak in Numerical Recipe (p 491, 3rd Ed.)
  """
  func = wrap_function(func,args)

  gold = GoldenRatio+1.0
  tiny = 1.e-20

  fa = func(ax)
  fb = func(bx)
 
  if fb > fa:   # switch roles of a and b so that the function goes downhill in the a -> b direction 
    t = ax; ax = bx; bx = t
    t = fa; fa = fb; fb = t
  
  cx = bx + gold*(bx-ax)
  cx = min(cx, xmax)
  cx = max(cx, xmin) 
  fc = func(cx)

  u_old = 1.e10
  if debug:
    print "x=%12.6f, f(x)=%16.6f" %(ax,fa)
    print "x=%12.6f, f(x)=%16.6f" %(bx,fb)
    print "x=%12.6f, f(x)=%16.6f" %(cx,fc)
  n_eval = 3

  if save_all:
    save=[(ax,fa),(bx,fb),(cx,fc)]

  while fb > fc: 
    r = (bx-ax)*(fb-fc)
    q = (bx-cx)*(fb-fa)
    u = bx-( (bx-cx)*q-(bx-ax)*r)/ (2.*copysign(max(abs(q-r),tiny),q-r))

    if cx > bx:      # set the limit of u
      if xmax > 1.e10: 
        ulim = bx + glimit*(cx-bx)
      else:
        ulim = xmax
    else:
      if xmin < -1.e10:
        ulim = bx + glimit*(cx-bx)
      else:  
        ulim = xmin

    if (bx-u)*(u-cx) > 0.0:   # u in [b,c]
      fu = func(u)
      n_eval += 1
      if save_all: save.append( (u,fu) )
      if debug: print "x=%12.6f, f(x)=%16.6f" %(u,fu)

      if fu < fc:  # minimum in [b,c]
        ax=bx; bx=u; fa=fb; fb=fu
        break 
      elif fu> fb:  # minimum in [a,u]
        cx=u; fc=fu
        break 

      u = cx + gold*(cx-bx)
      u = min(u, xmax)
      u = max(u, xmin) 
      fu = func(u)
      n_eval += 1
      if save_all: save.append( (u,fu) )
      if debug: print "x=%12.6f, f(x)=%16.6f" %(u,fu)

    elif (cx-u)*(u-ulim) > 0.0 : # parabolic fit gives u in [cx,ulim]
      fu = func(u)
      if fu < fc:
        (bx,cx,u)  = (cx,u,u+gold*(u-cx))
        (fb,fc,fu) = (fc,fu,func(u))
        n_eval += 1
        if save_all: save.append( (u,fu) )
        if debug: print "x=%12.6f, f(x)=%16.6f" %(u,fu)

    elif (u-ulim)*(ulim-cx) >= 0.0:
      u = ulim
      fu = func(u)
      if save_all:
        save.append( (u,fu) )
      n_eval += 1
    else:
      u = cx+gold*(cx-bx)
      fu = func(u)
      n_eval += 1
      if save_all: save.append( (u,fu) )
      if debug: print "x=%12.6f, f(x)=%16.6f" %(u,fu)

    (ax,bx,cx) = (bx,cx,u)
    (fa,fb,fc) = (fb,fc,fu)

    if ax==bx or bx==cx or u == u_old: 
      print "ERROR: search is trapped to the same value !!!"
      print "  -- a possible cause is that [xmin,xmax] does not bracket a minimum"
      sys.exit(1)

    u_old = u

  if save_all:
    return [ax,bx,cx],[fa,fb,fc],n_eval,save
  else:
    return [ax,bx,cx],[fa,fb,fc],n_eval      

def f_Min1D_BracketFit(x0,fx0,dx,func,args=(),fit_type='H',xmin=-1.1e10,xmax=1.1e10,fit_tol=0.1,glimit=10.0,debug=False):
  """
  This function first bracket the interval in which a minimum exists, and then determine the minimal point by a fitting
  """
  func = wrap_function(func,args)

  gold = 1.618
  tiny = 1.e-20

  ax = x0
  fa = fx0
  bx = x0+dx
  fb = func(bx)
 
  if fb > fa:   # switch roles of a and b so that the function goes downhill in the a -> b direction 
    t = ax; ax = bx; bx = t
    t = fa; fa = fb; fb = t

  cx = bx + gold*(bx-ax)
  cx = min(cx, xmax)
  cx = max(cx, xmin) 
  fc = func(cx)
  n_eval=2
  
  xdata=[ax,bx,cx]
  fdata=[fa,fb,fc]

  u_old = 1.e10
  if debug:
    save=[(ax,fa),(bx,fb),(cx,fc)]

  while fb > fc: 
    r = (bx-ax)*(fb-fc)
    q = (bx-cx)*(fb-fa)
    u = bx-( (bx-cx)*q-(bx-ax)*r)/ (2.*copysign(max(abs(q-r),tiny),q-r))

    if cx > bx:      # set the limit of u
      if xmax > 1.e10: 
        ulim = bx + glimit*(cx-bx)
      else:
        ulim = xmax
    else:
      if xmin < -1.e10:
        ulim = bx + glimit*(cx-bx)
      else:  
        ulim = xmin

    if (bx-u)*(u-cx) > 0.0:   # u in [b,c]
      fu = func(u)

      if fu < fc or fu > fb: 
        xdata.append(u)
        fdata.append(fu)
        n_eval += 1
        if debug: save.append( (u,fu) )
        if fu < fc:  # minimum in [b,c]
          ax=bx; bx=u; fa=fb; fb=fu; break 
        elif fu> fb:  # minimum in [a,u]
          cx=u; fc=fu; break 

      u = cx + gold*(cx-bx)
      u = min(u, xmax)
      u = max(u, xmin) 
      fu = func(u)
    elif (cx-u)*(u-ulim) > 0.0 : # parabolic fit gives u in [cx,ulim]
      fu = func(u)
      if fu < fc:
        (bx,cx,u)  = (cx,u,u+gold*(u-cx))
        (fb,fc,fu) = (fc,fu,func(u))
    elif (u-ulim)*(ulim-cx) >= 0.0:
      u = ulim
      fu = func(u)
    else:
      u = cx+gold*(cx-bx)
      fu = func(u)

    xdata.append(u)
    fdata.append(fu)
    n_eval += 1
    if debug: 
      save.append( (u,fu) )

    (ax,bx,cx) = (bx,cx,u)
    (fa,fb,fc) = (fb,fc,fu)

    if ax==bx or bx==cx or u == u_old: 
      print "ERROR: search is trapped to the same value !!!"
      print "  -- a possible cause is that [xmin,xmax] does not bracket a minimum"
      sys.exit(1)

    u_old = u
  # end the bracketing loop  

  if debug:
    print "\n\tBracketed intervals:"
    print "\tax=%12.6f, fa=%12.6f"%(ax,fa)
    print "\tbx=%12.6f, fb=%12.6f"%(bx,fb)
    print "\tcx=%12.6f, fc=%12.6f\n"%(cx,fc)

  # a minimum is bracketed, xmin and fmin are determined by a fitting 
  fit_par,fit_err,xy0 = f_EOS_Fitting(xdata,fdata,fit_type)

  # the followings make some "security" checks: 
  # if the fitted x_min is not in the bracketed interval, 
  # or if f(x_min) is greater than the minimal f(x) in the scanned range 
  # of if the fitting error is too big, 
  #  the fitted value is discarded. Instead the value of x that gives minimal f(x) in 
  #  the scanned range is returned 
  
  xm = fit_par[0]
  fm = fit_par[1]

  if debug:
    print "\n\tData to be fitted:"
    print "\t%6s %12s %16s"%("i","x","f(x)") 
    for i in range(len(xdata)):  print "%6d %12.6f %16.6f" %(i,xdata[i],fdata[i])
    print "\n\txm(fit)= %12.6f, fm(fit)=%12.6f\n"%(xm,fm)

  x_start =  min(xdata)
  x_end   =  max(xdata)

  if xm < x_start or xm > x_end or fm > min(fdata) or fit_err > fit_tol :
    status = 1
    if fm > min(fdata):
      print "WARNING: fmin from fitting is bigger than one of the fitted values" 

    if xm < x_start or xm > x_end:
      print "WARNING: x_min found from fitting falls out of the bracketed interval (%12.6f,%12.6f)" %(x_start,x_end)

    if fit_err > fit_tol:
      print "WARNING: fitting error=%12.2e is bigger than the tolerance (%12.2e)" %(fit_err,fit_tol)

  else: 
    status = 0
    xdata.append(xm) 
    f_xm =  func(xm)
    fdata.append(f_xm)
    if debug: print "\tf(xm_fit) = %12.6f"%(f_xm) 

  k_min = array(fdata).argmin()
  x_min = xdata[k_min]
  f_min = fdata[k_min]

  if status ==0 and x_min != xm: 
    print "WARNING: x_min != xm(fit)  !!!"

  if debug: print"\txm=%12.6f, f(xm)=%12.6f)\n"%(x_min,f_min)

  if debug:
    return x_min,f_min,x_start,x_end,status,fit_par,fit_err,save
  else:
    return x_min,f_min,x_start,x_end,status      

def f_Min1D_Brent(func,xbrak,fbrak=[],args=(),xtol=1.e-4,ftol=1.e-4,itmax=100,debug=False,save_all=False):
  """
  This is a modified implemetantion of Brent's method as described in Numerical Recipe (3rd Ed. pp 496)
  """
  zeps=1.e-10
  cgold=0.381966

  n_eval = 0

  func = wrap_function(func,args)

  bx = xbrak[1] 
  if xbrak[0] < xbrak[2]:
    ax = xbrak[0]
    cx = xbrak[2]
  else: 
    ax = xbrak[2]
    cx = xbrak[0]

  if not fbrak:
    fa=func(ax)
    fb=func(bx)
    fc=func(cx)
    n_eval += 3
  else: 
    fb = fbrak[1]
    if xbrak[0] < xbrak[2]:
      fa = fbrak[0]
      fc = fbrak[2]
    else:
      fa = fbrak[2]
      fc = fbrak[0]

  a = ax
  b = cx 
  x = bx; fx = fb
  if fa < fc:
    w = ax; fw = fa
    v = cx; fv = fc
  else:
    w = cx; fw= fc
    v = ax; fv = fa

  u = x;  fu = fx
  
  d = 0.0
  if x>= 0.5*(a+b):
    e = a-x
  else:
    e = b-x

  if save_all: save=[(ax,fa),(bx,fb),(cx,fc)]
  
  for iter in range(itmax): 
    xm=0.5*(a+b)
    tol1=xtol*abs(x)+zeps
    tol2=2.0*tol1

    u_old = u
    fu_old = fu

    if abs(x-xm) <=  (tol2-0.5*(b-a)):  break

    if  abs(e) > tol1: 
      r=(x-w)*(fx-fv)
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=2.0*(q-r)
      if q > 0.0: p = -p
      q=abs(q)
      etemp=e
      e=d
      if abs(p) >= abs(0.5*q*etemp) or p <= q*(a-x) or p >= q*(b-x):
        if x>= xm: 
          e = a-x
        else:
          e = b-x
        d = cgold*e

      else:
        d=p/q; u=x+d
        if u-a < tol2 or b-u < tol2: 
          d=copysign(tol1,xm-x)

    else:
      if x>= xm : e = a-x 
      else: e = b-x
      d = cgold*e

    if abs(d) >= tol1:
      u = x+d
    else:
      u = x+copysign(tol1,d)

    fu=func(u)
    n_eval += 1
    if debug: print "iter=%d, x=%12.6f f(x)=%16.6f" %(iter,u,fu)
    if save_all: save.append((u,fu))
    
    if fu <= fx: 
      if u >= x: a=x
      else: b=x
      (v,w,x)=(w,x,u)
      (fv,fw,fx) = (fw,fx,fu)
    else:
      if u < x: a=u
      else: b=u
      if fu <= fw or w == x:
        v=w; w=u; fv=fw; fw=fu
      elif fu <= fv or v == x or v == w:
        v=u; fv=fu

    if abs(fu_old-fu) < ftol or abs(u_old-u) < xtol: break
      
  # end of loop over iter 
  if iter >= itmax:
    print "WARNING: Convergence fail to reach in f_Min1D_Brent "

  if save_all:
    return (x,fx,n_eval,save)
  else: 
    return (x,fx,n_eval)  

#
# This defines the 1-D function used in the line minimization of N-dimensional optimzation problem   
#
def f_LinMin_F1D(x,func,pvec,nvec,args_func):
  """
  This defines the 1-D function used in the line minimization of N-dimensional 
  optimzation problem
  """
  func = wrap_function(func,args_func)

  px = pvec.copy()
  px += x * nvec
  return func(px.tolist())


   
def f_LinMin(func_nd,pvec,nvec,xtol,ftol,args_func_nd,bounds=(-1.e20,1.e20),save=[],save_all=False):
  """
  This function perform a line-minimization F( pvec + x*nvec) 
  """
  ax=0.0; bx=1.0
  n_eval=0
  args_f1d = (func_nd,pvec,nvec,args_func_nd)
  f1d = f_LinMin_F1D 
  ndim = pvec.shape[0]

  xmin,xmax = bounds

  if save_all :
    xbrak,fbrak,n_eval_brak,save_1d = f_Min1D_Bracket(f1d,ax,bx,xmin=xmin,xmax=xmax,args=args_f1d,save_all=save_all)

    for i in range( len(save_1d) ):
      x = save_1d[i][0]
      fx = save_1d[i][1]
      ptmp = pvec + x* nvec
      save.append( (ptmp.tolist(),fx) )
     
    xmin,fmin,n_eval_min,save_1d = f_Min1D_Brent(f1d, xbrak, fbrak=fbrak,args=args_f1d,xtol=xtol,ftol=ftol,save_all=save_all)
    for i in range( len(save_1d) ):
      x = save_1d[i][0]
      fx = save_1d[i][1]
      ptmp = pvec + x* nvec
      save.append( (ptmp.tolist(),fx) )

  else:
    xbrak,fbrak,n_eval_brak = f_Min1D_Bracket(f1d,ax,bx,xmin=xmin,xmax=xmax,args=args_f1d)
    xmin,fmin,n_eval_min = f_Min1D_Brent(f1d, xbrak, fbrak=fbrak, args=args_f1d, xtol=xtol,ftol=ftol)

  n_eval += n_eval_brak + n_eval_min
  nvec *= xmin 
  pvec += nvec

  return fmin,n_eval 

def f_Minimize_Powell_NR(func,x0,args=(),xtol=1.e-4,ftol=1.e-5,itmax=200,debug=False,save_all=False):
  """
  This function optimizes a n-d function by using Modified Powell's method as presented in Numerical Recipe
  """

  func = wrap_function(func,args)
  tiny=1.e-10
  p=array(x0)
  ndim = p.shape[0]
  save = []
  n_eval = 0
  print "ndim=%d" %(ndim)
  
  fmin=func(p.tolist())
  n_eval += 1
  if save_all:
    save.append( (p.tolist(),fmin) )

  xi = zeros( (ndim,ndim) )
  for i in range(ndim):
    xi[i,i]=1.0    
    
  iter = 0
  pt = p.copy()      # save the initial point 
  while iter < itmax:
    iter += 1
    fp = fmin
    p_old = p.copy()
    ibig = 0
    decl = 0.

    for i in range(ndim):
      xit = xi[i].copy()
      fptt = fmin 

      fmin,n_eval_line = f_LinMin(func, p, xit, xtol,ftol, args,save_all=save_all,save=save) 
      n_eval += n_eval_line
      if debug: print "\tcalling func # %d" %(n_eval_line)

      if fptt-fmin > decl: 
        decl = fptt-fmin 
        ibig = i 
    # end loop over i
    
    ferr = abs(fp-fmin)
    xerr = max( abs(p_old-p) )
    if debug:
      print "i=%6d, xerr= %12.6f, ferr=%12.6f" %(iter,xerr,ferr)
    if ferr < ftol or xerr < xtol: break 

    ptt = 2.*p - pt
    xit = p - pt
    pt  = p.copy()
    
    fptt = func(ptt.tolist())
    n_eval += 1
    if save_all:
      save.append( (ptt.tolist(),fptt) )

    if fptt >= fp: continue 
    t = 2.0*(fp-2.0*fmin+fptt)*(fp-fmin-decl)**2-decl*(fp-fptt)**2
    if t >= 0.0 : continue 
    
    fmin,n_eval_line = f_LinMin(func, p, xit, xtol,ftol, args,save_all=save_all,save=save)
    if debug: print "\tcalling func # %d in line minimization" %(n_eval_line)
    n_eval += n_eval_line

    xi[ibig] = xi[ndim-1].copy()
    xi[ndim-1] = xit.copy()
    
  # end of loop over iter 
  if iter >= itmax:
    print "ERROR in f_Powerll: convergence not reached after %d iterations"%(itmax)
  
  if save_all:
    return p.tolist(),fmin,save
  else:
    return p.tolist(),fmin

def f_Minimize_Powell(func,x0,args=(),xtol=1.e-4,ftol=1.e-5,itmax=200,debug=False,save_all=False,restart=False):
  """
  This function optimizes a n-d function by using Powell's method 
  """

  tiny=1.e-10
  ndim = len(x0)
  save = []
  n_eval = 0

  x = x0[:]; fx = 0.0
  if restart: x,fx=f_Read_Restart()

  # set up the initial directional matrix 
  xi = zeros( (ndim,ndim) )
  for i in range(ndim): xi[i,i]=1.0    
    
  iter = 0
  p=array(x); fp = fx
  while iter < itmax:
    iter += 1
    fp0 = fp
    p0 = p.copy()

    # write out results of the current cycyle to be used for restarting 
    ofile = open(opt_restart_file,'a')
    for i in range(ndim): ofile.write(" %24.12e " %(p[i]))
    ofile.write(" %24.12e \n" %(fp))
    ofile.close()

    for i in range(ndim):
      xit = xi[i].copy()
      fp,n_eval_line = f_LinMin(func,p,xit,xtol,ftol,args,save_all=save_all,save=save) 

      f_Write_Restart(p.tolist(),fp)

      n_eval += n_eval_line
      if debug: print "\tcalling func # %d" %(n_eval_line)
    # end loop over i

    if ndim == 1: break 

    ferr = abs(fp-fp0)
    xerr = max(abs(p-p0))
    if debug:
      print "  iter=%d,x_err=%12.6e, f_err=%12.6e" %(iter,xerr,ferr)
    if ferr < ftol or xerr < xtol : break
   
    for i in range(ndim-1): xi[i]=xi[i+1].copy()
    xi[ndim-1]= p - p0

    xit = xi[ndim-1].copy()
    p= p0.copy() 
    fp,n_eval_line = f_LinMin(func,p,xit,xtol,ftol, args,save_all=save_all,save=save)

    f_Write_Restart(p.tolist(),fp)

    if debug: print "\tcalling func # %d in line minimization" %(n_eval_line)
    n_eval += n_eval_line
    
  # end of loop over iter 
  if iter >= itmax:
    print "ERROR in f_Powerll: convergence not reached after %d iterations"%(itmax)
  
  if save_all:
    return p.tolist(),fp,save
  else:
    return p.tolist(),fp

def f_Minimize_TaxiCab(func,x0,\
                       args=(),\
                       constraint = None, \
                       bounds=[],\
                       xtol=1.e-4,\
                       ftol=1.e-5,\
                       itmax=200,\
                       debug=False,\
                       save_all=False,\
                       restart=False):
  """
  This function optimizes a n-dimensional function by using the Taxi Cab method 
  """
  tiny=1.e-10
  ndim = len(x0)
  save=[]

  if constraint is None: constraint = f_List_New(1,ndim) 

  x = x0[:]; fx = 0.0
  # enable restart mode
  if restart: x,fx=f_Read_Restart()

  # set up the unitary n X n directional matrix 
  xi = zeros( (ndim,ndim) )
  for i in range(ndim): 
    xi[i,i]=1.0    
    
  p=array(x)
  fmin = fx
  iter = 0
  while iter < itmax:
    iter += 1
    p_old = p.copy()
    fp = fmin
    ibig = 0
    decl = 0.
   
    for i in range(ndim):
      if constraint[i] == 0: continue 

      xit = xi[i].copy()
      if len(bounds) > 0 :
        bounds_i = bounds[i]
      else:
        bounds_i = (-1.e20,1.e20)

      fmin,n_eval_line = f_LinMin(func,p,xit,xtol,ftol,args,bounds = bounds_i, save_all=save_all,save=save) 
      if debug: print "\tcalling func # %d in line minimization" %(n_eval_line)

      f_Write_Restart(p.tolist(),fmin)

    # end loop over i

    if ndim == 1: break 

    ferr = abs(fp-fmin)
    xerr = max(abs(p-p_old))
    if debug:
      print "  iter=%d,x_err=%12.6f, f_err=%12.6f" %(iter,xerr,ferr)
    if ferr < ftol or xerr < xtol : break 
    
  # end of loop over iter 
  if iter >= itmax:
    print "ERROR in f_Minimize_TaxiCab: convergence not reached after %d iterations"%(itmax)
  
  if save_all:
    return p.tolist(),fmin,save
  else:
    return p.tolist(),fmin


def f_Minimize_ScanFit(func,x0,xranges,\
                       constraint=None, \
                       args=(),\
                       xtol=1.e-4,\
                       ftol=1.e-6,\
                       fit_tol=1.e-1,\
                       fit_types=[],\
                       itmax=10,\
                       debug=False,\
                       save_all=False,\
                       restart=False):
  """
  This function optimizes a n-dimensional function by using the scan-and-fit approach
  i.e. each variable is scanned in a given range and the function is fitted by 
  a function as assigned by fit_types 
  """
  ndim = len(x0)
  func = wrap_function(func,args)

  if constraint is None: constraint=f_List_New(1,ndim) 

  save=[]

  x = x0[:]; fx = 0.0
  if restart: x, fx = f_Read_Restart()
  
  iter = 0
  while iter < itmax: 
    iter += 1
    x_old = x[:]
    fx_old = fx 

    # write out results of the current cycyle to be used for restarting 
    f_Write_Restart(x,fx)

    for i in range(ndim):
      if constraint[i] == 0: continue 

      ## the scan range is designated with respect to the current "mimimal" point 
      xi_start = x[i] + xranges[i][0]
      xi_end = x[i] + xranges[i][1]
      if xi_start > xi_end: 
        t = xi_start; xi_start = xi_end; xi_end = t

      xi_nsteps= xranges[i][2]    
      xi_step = (xi_end - xi_start)/(xi_nsteps-1)
 
      fxi_data =[]
      xi_data =[]
      if debug:
        print "\n  Scan %d-th variable:\n\t%12s %16s"%(i,"x"+str(i),"f(x)")
      for k in range(xi_nsteps): 
        xi = xi_start + k*xi_step
        x[i] = xi
        fx = func(x)
        xi_data.append(xi)
        fxi_data.append(fx)
        if save_all: save.append( (x,fx) )
        if debug:
          print "\t%12.6f %16.6f" %(xi,fx)
      
      if len(fit_types) == 0: fit_type = 'H' 
      else: fit_type = fit_types[i]
      fit_par,fit_err = f_EOS_Fitting(xi_data,fxi_data,fit_type)

      if debug:
        print "  Fitting results for %d-th variable in %d-th cycle :" %(i,iter)
        print "\tfitting method: " + fit_type
        print "\tfitted parameters:",fit_par
        print "\tfitting error = %12.2e:" %(fit_err)

      # check whether something abnormal is occurring: ideally, the minimum should fall within the scanned range. 
      # but occassioally that is not the case. That might be OK if the function decreases at the fitted minimum  
      # the followings make some "security" checks: first whether the fitted minimum is in the scan range, and second whether 
      # the function at the fitted minimum is smaller than any of the scanned values. 
      # If the latter is not the case, then the fitted value is discarded and the scanned value that 
      # gives the lowest function value is used instead   
      # if the fitting error is too big, the fitted value is also discarded  
      xi_min = fit_par[0]
      fxi_min = fit_par[1] 
      x[i] = xi_min 
      fx = fxi_min
     
      xi_min_out_of_range = False
      if xi_min < xi_start or xi_min > xi_end:
        xi_min_out_of_range = True
      
      fxi_min_increasing = False
      if fxi_min > array(fxi_data).min():
        fxi_min_increasing = True

      if fxi_min_increasing or fit_err > fit_tol:
        k_min = array(fxi_data).argmin()
        x[i] = xi_data[k_min]
        fx = fxi_data[k_min]
        
      if fxi_min_increasing:
        print "WARNING: fmin from fitting is bigger than one of the fitted value when scanning %d-th variable" %(i)
        print "  -- fitted parameters are discarded!!!"
    
      if xi_min_out_of_range: 
        print "WARNING: x_min found from scan&fit %d-th variable falls out of the scanned range (%12.6f,%12.6f)" %(i,xi_start,xi_end)

      if fit_err > fit_tol:
        print "WARNING: fitting error is bigger than the tolerance (%12.6f)" %(fit_tol)
        print " -- fitted parameters are discarded!!!"
        
    ## end of loop over i   
    # check convergence 


    ferr = abs( fx - fx_old )
    xerr = 0.0
    for i in range(ndim):
      xi_err = abs(x[i]-x_old[i])
      if xi_err > xerr: xerr = xi_err

    if debug:
      print "  iter=%d, x_err=%12.6f, f_err=%12.6f" %(iter,xerr,ferr)

    if xerr < xtol or ferr < ftol: break 
  # end of loop over iter 
  if iter >= itmax:
    print "WARNING: Fail to converge by the Scan&Fit method"

  fx = func(x)
  if save_all: save.append( (x,fx) )
    
  if save_all:
    return x,fx,save
  else:
    return x,fx 

#
# This defines the 1-D function used in the Bracket&Fit minimization for N-dimensional optimzation problem   
#
def f_Func1D(x_i,func_nd,x,i,args_func):
  """
  This defines the 1-D function used in the line minimization of N-dimensional 
  optimzation problem
  """
  func = wrap_function(func_nd,args_func)

  xt=x[:]
  xt[i] = x_i
  return func(xt)

def f_Minimize_BracketFit(
                       func, x0,\
                       constraint=None,\
                       dxs=[],\
                       args=(),\
                       xtol=1.e-4,\
                       ftol=1.e-6,\
                       fit_tol=1.e-2,\
                       fit_types=[],\
                       itmax=50,\
                       debug=False,\
                       restart=False):
  """
  This function optimizes a n-dimensional function by using the bracket-and-fit approach
  i.e. an approximate minimum of each variable is first bracketed and then the function is fitted by 
  a parabolic function 
  """
  func_nd = wrap_function(func,args)
  func_1d = f_Func1D

  ndim = len(x0)

  # constraint is an integer array to control whether a particular variable should be fixed 
  # (if set 0) 
  if constraint is None : constraint = f_List_New(1,ndim)

  # dxs(1:ndim) -- the size of initial step 
  if len(dxs) == 0 : dxs=f_List_New(2.0,ndim)

  # the type of functions used for fitting, the default is a harmonic function 
  if len(fit_types) == 0: fit_types=f_List_New('H',ndim)

  if debug: save_all = True

  if restart: 
    x, fx = f_Read_Restart()
  else: 
    x = x0[:]; fx = func_nd(x)

  iter = 0
  while iter < itmax: 
    iter += 1
    x_old = x[:]
    fx_old = fx 

    # write out results of the current cycyle to be used for restarting 
    f_Write_Restart(x,fx)

    if debug:
      msg = "\n\nicycle=%6d "%(iter) 
      for i in range(ndim): msg = msg+"%12.6f "%(x[i])
      msg = msg + "%12.6f "%(fx)
      print msg

    for i in range(ndim):
      if constraint[i] == 0:  
        print "  Constrained optimization is required: %d-th variable is skipped " %(i)
        continue 
 
      print "  Start bracket&fit over %d-th variable in %d-th cycle :" %(i,iter)
      args_1d = (func,x,i,args)   # arguments for 1-D function 
      x0_i = x[i]
      fx0_i = fx
      dx_i = dxs[i]
      fit_type_i= fit_types[i]
      out = f_Min1D_BracketFit(x0_i,fx0_i,dx_i,func_1d,args=args_1d,fit_type=fit_type_i,fit_tol=fit_tol,debug=debug)

      x[i] = out[0]
      fx   = out[1]
      status = out[4]

      # if the bracket-fitting is successful, as the iteraion moves towards the minimum, 
      # the inital guess for the braketing should also 
      # decrease so that the harmonic fitting can be come more accurate 
      # dxi_start and dx_end measures the distance between the current xi_min and the bracketing bounds
      
      if status == 0: 
        dxi_start = out[2]-x[i]
        dxi_end = out[3]-x[i]
        if abs(dxi_start) > abs(dxi_end): 
          dxs[i] = dxi_end/2
        else:
          dxs[i] = dxi_start/2
      
      if debug:
        print "  Fitting results:" 
        print "\tFitting method: " + fit_type_i
        print "\tx_min=%12.6f,f_min =%20.6f " %(x[i],fx)
        print "\tthe bracketed interval=(%12.6f,%12.6f)"%(out[2],out[3])  
        print "\tFitting parameters:",out[5]
        print "\tFitting error = %12.2e:" %(out[6])
        print "\tTotal number of function evaluations:%6d\n" %(len(out[7]))

    ## end of loop over i   
    # check convergence 

    ferr = abs( fx - fx_old )
    xerr = 0.0
    for i in range(ndim):
      xi_err = abs(x[i]-x_old[i])
      if xi_err > xerr: xerr = xi_err

    if debug:
      print "  iter=%d, x_err=%12.6f, f_err=%12.6f" %(iter,xerr,ferr)

    if (status == 0 ) and (xerr < xtol or ferr < ftol): break 

  # end of loop over iter 
  if iter >= itmax: print "WARNING: Fail to converge by the Scan&Fit method"

  # perform a final function evaluation: this is necessary if this function is 
  # used for wien2k lattice optimzation so that the current directory contains lattice 
  # structure that correspond to the optimal lattice structures 
  fx = func_nd(x)

  return x,fx 
   
