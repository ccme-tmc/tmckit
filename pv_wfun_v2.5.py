#!/usr/bin/env python

# v1 is written by HJ, v2 is modified by ZHC

import sys,os,shutil,scipy,numpy,math
import matplotlib.pyplot as plt
from scipy import interpolate
from struct_utils import *
from constants import *
#from vasp_utils_jh import vasp_read_pot
from io_utils import *
#from wfun_utils import *
#import matplotlib.pyplot as plt

def vasp_read_pot(potfile="LOCPOT",zshift=0.0,bulk_flag=0,direct_flag='0'):
  """
  Read the electron density or potential from the file 
  """
  ifile = open(potfile,'r')
  ifile.readline()

  latt_scale = float(ifile.readline())  # scale for the lattice vectors

  # read the lattice vectors 
  latt_vec=[]
  for i in range(3):
    line_s = ifile.readline().split()
    x = float(line_s[0])
    y = float(line_s[1])
    z = float(line_s[2])
    latt_vec.append([x,y,z])
  #print latt_vec 
  mat=numpy.array(latt_vec)
  if direct_flag=='0':
	  if mat[2][2]>mat[0][0]:
		if mat[2][2]>mat[1][1]:
			zdirect=2
		else:
			zdirect=1
	  else:
		if mat[0][0]>mat[1][1]:
			zdirect=0
		else:
			zdirect=1
  else:
	  if direct_flag=='x':
	  	zdirect=0
	  elif direct_flag=='y':
	  	zdirect=1
	  else:
	  	zdirect=2
#  print mat
  trans_mat=numpy.transpose(mat) 
 # print trans_mat
  inv_mat=numpy.linalg.inv(trans_mat)
 # print inv_mat
  # read species

  specs = ifile.readline().split()

  nat_sp =[]
  line_s = ifile.readline().split()
  nat = 0
  for i in range(len(specs)):
    nat_sp.append(int(line_s[i]))
    nat += nat_sp[i]
  
  cord_type=ifile.readline()
  #add by ZHC
  cord_type=cord_type.lstrip()
  zlist=[]	
  zmin=1000000.0
  zmax=0.0
  for i in range(nat):
	line = ifile.readline().split()
	atvec=[]
	for j in range(3):
		atvec.append(float(line[j]))
	atvec=numpy.array(atvec)
#	print atvec
#	zpos=float(line[2])
	if cord_type[0]=='D' or cord_type[0]=='d':                          ### all convert to Cartizian to find the center of the slab
		#print "NOTE: The coordination is DIRECT!!!"	
		for j in range(3):
			if j==zdirect:
				atvec[j]=atvec[j]+zshift
			atvec[j]=atvec[j]-math.floor(atvec[j])		
		cart_vec=numpy.dot(trans_mat,atvec)
		if cart_vec[zdirect]<zmin:
			zmin=cart_vec[zdirect]
		if cart_vec[zdirect]>zmax:
			zmax=cart_vec[zdirect]	
	else:
		#print "NOTE: The coordination is CARTISIAN!!!"
		dir_vec=numpy.dot(inv_mat,atvec)
		#print dir_vec
		for j in range(3):
			if j==zdirect:
				dir_vec[j]=dir_vec[j]+zshift
			dir_vec[j]=dir_vec[j]-math.floor(dir_vec[j])
		#print dir_vec
		
		cart_vec=numpy.dot(trans_mat,dir_vec)
		#print cart_vec
		if cart_vec[zdirect]<zmin:
			zmin=cart_vec[zdirect]
		if cart_vec[zdirect]>zmax:
			zmax=cart_vec[zdirect]	
 	
  center_pos=(zmax+zmin)/2.0

#  print "zmax,zmin"
#  print zmax,zmin
#  print "CENTER:"
#  print center_pos
  dis_vac_up=mat[zdirect][zdirect]-zmax
  dis_vac_down=zmin-0.0
  thick_vac=abs(dis_vac_up+dis_vac_down)
  
  if zmax-zmin>0.8*mat[zdirect][zdirect] and mat[zdirect][zdirect]>12.0 and bulk_flag==0:
	frac=(mat[zdirect][zdirect]-zmax)/float((mat[zdirect][zdirect]))
	print "PLEASE SET -z tag to ensure the slab does not cross the boundary!!!"
	sys.exit(0)	
  vac_center=thick_vac/2.0+zmax
  if vac_center>mat[zdirect][zdirect]+zshift*mat[zdirect][zdirect]:
	vac_center=vac_center-mat[zdirect][zdirect]
#  print "vac_center" 
#  print vac_center
 
  ifile.readline()

  line_s = ifile.readline().split()
  nx = int(line_s[0])
  ny = int(line_s[1])
  nz = int(line_s[2]) 
  ndata = nx*ny*nz
  nlines = ndata/5
  if ndata%5 != 0: nlines += 1
  data = []
  for il in range(nlines):
    line_s = ifile.readline().split()
    for i in range(len(line_s)):
      data.append(float(line_s[i]))
  ifile.close()
  return data,(nx,ny,nz),latt_vec,center_pos,vac_center,zmax,zmin,zdirect,thick_vac


#vasp_read_pot()

def f_Help_Info():
  myname = sys.argv[0]
  print "\n"+myname + ": a Python script to get work function using VAPS output. Mainly for bulk correction calculations. \n"
  print "   Usage: " + os.path.basename(myname) + " [options]"\
  + """
  Options:
    -h                # display this help information
    -f <file >        # prefix for input and output
    -d <x/y/z>        # the normal direction of the surface. Default is determined automatically.
    -n <int>          # set the # of point for average WITHOUT automatic procedure. Use for the case that the average range is not correctly found.
    -i <int>          # set the factor of cubic spline interpolation. The default value is 10, so that the total number of point will be 10*number of original data. If i=0, the interpolation will be switched off.
    -b <1/0>          # set the type slab(0), bulk(1). Default: bulk or slab will be determined automatically. Or you can use this tag to force it as bulk.
    -t <data_type>:   # VH/RHO
    -tol <float>:     # the tolerance for automatic determination of periodic range for averaging. The default value is 0.2.
    -z <float>:       # shift of the fractional coordination in z direction to ensure the slab will NOT span the boundary. The default value is 0.0.
    -s <int>          # 1 for save the electric potential curve figure as png. 0 for show the figure to the screen. The default is 0.
    --debug           # set debug mode
  Examples1:
	./pv_wfun_v2.py -f LOCPOT
  """
  sys.exit(0)

def findpeak(ax,ay):                                           #find the peaks of electric potential
	if len(ax)!=len(ay):
		print "FIND PEAK ERROR: length of vector not equal"
	peakx=[]
	peaky=[]
	peakindex=[]
	for i in range(len(ay)-1):
		if (ay[i]-ay[i-1])*(ay[i+1]-ay[i])<=0.0:
			peakindex.append(i)	
			peakx.append(ax[i])
			peaky.append(ay[i])
	last=len(ay)-1
	if (ay[last]-ay[last-1])*(ay[0]-ay[last])<=0.0:
		peakx.append(ax[last])
		peaky.append(ay[last])
		peakindex.append(last)
	if len(peakx)<=3:
		print "\n WARNING!!!!!!!!!!!!\n The # of peaks is rather small!!!!!The averaging range may be NOT correct!!!!\n"
	return peakx,peaky,peakindex

# default input and output format 

def_dnorm = '0'
#def_intflag = 1
def_file = ''
def_debug = False
def_np = 10 
def_bulk=1
def_zshift=0.0
def_tolerance=0.2
def_save=0


if f_Getopt('-h',0,False): f_Help_Info()

dnorm = f_Getopt('-d',1,'0' )
np = f_Getopt('-i',1,10) 
bulk = f_Getopt('-b',1,0) 
num_set=f_Getopt('-n',1,0)
file  = f_Getopt('-f',1,'LOCPOT' ) 
#np    = f_Getopt('-n',1,10   ) 
type  = f_Getopt('-t',1,'VH') 
zshift  = f_Getopt('-z',1,0.0) 
tolerance  = f_Getopt('-tol',1,0.2) 
save = f_Getopt('-s',1,0) 



vh,nxyz,latt_vec,center_pos,vac_center,zmax,zmin,zdirect,thick_vac = vasp_read_pot(file,zshift,bulk,dnorm)

if bulk==0 and thick_vac<5.0:
	bulk=1



directs=['x','y','z']
if dnorm=='0':
	dnorm=directs[zdirect]

#print center_pos
(nx,ny,nz) = nxyz
#print " nx,ny,nz=",nx,ny,nz

if dnorm == 'x':
  nd = nx
  vd = latt_vec[0]
elif dnorm == 'y':
  nd = ny
#  vd = latt_vec[1+shift]
  vd = latt_vec[1]
elif dnorm == 'z':
  nd = nz
  vd = latt_vec[2]
elif dnorm == '0':
  nd = 1
  vd = [0.0,0.0,0.0]
else:
  print "ERROR: unsupported option for dnorm= ",dnorm

nsum = nx*ny*nz/nd


vd0 = []
for i in range(nd): vd0.append(0.0)

ip=0
vsum=0.0
for iz in range(nz):
  for iy in range(ny):
    for ix in range(nx):
      if dnorm == 'x': vd0[ix] += vh[ip]
      elif dnorm == 'y': vd0[iy] += vh[ip]
      elif dnorm == 'z': vd0[iz] += vh[ip]
      else: vd0[0] += vh[ip]
      ip += 1

for i in range(nd): vd0[i] = (vd0[i]/nsum) 
print 'Number of original data:',nd

if nd == 1:
  print "Averaged electrostatic potetial= %12.4f"%(vd0[0])
  sys.exit(0)  

ld = sqrt(vd[0]**2+vd[1]**2+vd[2]**2)
vd0.append(vd0[0])
nd += 1

hd = ld/(nd-1)
xd0 = []
for i in range(nd): xd0.append(i*hd)

fname = file+'-'+type+'-wfun'

print "Original data after xy averaging is written to ",fname+"0.dat"

ofile = open( fname+"0.dat",'w') 
xd_moved=[]
for i in range(len(xd0)):
#  	xd_moved.append(xd0[i]-ld/2.0)
  	xd_moved.append(xd0[i]+zshift*vd[zdirect])
	ofile.write('%12.6f %12.6f \n' %(xd0[i]-ld/2.0,vd0[i]))
ofile.close()


if bulk==1:     
	sum_bulk=0.0                                             ### bulk==1 case
	for i in range(len(xd0)):
		sum_bulk=sum_bulk+vd0[i]
	sum_bulk=sum_bulk/float(len(xd0))
	print "TREAT AS BULK!!!"
	
	print "The averaged electric potential is:"
	print sum_bulk
  	sys.exit(0)  
        

#print xd_moved

#print "V0(z=0  )= %12.4f "%vd0[0]  

if np > 0:                                                   ### interpolation in cubic spline
  print 'Factor of interpolated data:',np  
  np = (nd-1)*np+1   
  f_intp=interpolate.interp1d(xd_moved,vd0,kind='cubic')
  xd_intp=numpy.linspace(xd_moved[0],xd_moved[-1],num=np)
  y_intp=f_intp(xd_intp)
  peakx,peaky,peakindex=findpeak(xd_intp,y_intp)                   ### find the peak
  cent_center=100000.0
  cent_index=-999999
  for i in range(len(xd_intp)):
	dis=abs(xd_intp[i]-center_pos)
	if dis<cent_center:
		cent_center=dis
		cent_index=i
  
  dis_center=100000.0
  dis_index=-999999
  for i in range(len(peakx)):                                       ### find the peak nearest to the center of slab
	dis=abs(peakx[i]-center_pos)
	if dis<dis_center:
		dis_center=dis
		dis_index=i
#  print "DIS"
#  print dis_center,dis_index,peakx[dis_index]
  
  
  min_dis=9999999.0
  min_index=-9999999
  #tolerance=0.2
  for i in range(len(peakx)):                                      ### find another peak which has the same value of the peak nearest to the center
	if abs(peaky[i]-peaky[dis_index])<tolerance:
		if abs(peakx[i]-peakx[dis_index])<min_dis and i!=dis_index and abs(y_intp[peakindex[i]-int(np/50)]-y_intp[peakindex[dis_index]-int(np/50)])<tolerance:
			min_dis=abs(peakx[i]-peakx[dis_index])
			min_index=i
 # print min_dis,min_index
 # print "!!!!!Q"
 # print min_index,dis_index
  find_peak_flag=1
  if min_index==-9999999:
  	min_dis=9999999.0
	min_index=-9999999
	print "\nWARNING!!!!!!!!!!!!!!!!!!!!!!!!!\nCan not find the range for averaging automatically!!!\n TRY TO INCREASE THE TOLERANCE AND FIND IT AGAIN!!!\nWARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
	tolerance=2*tolerance
  	for i in range(len(peakx)):                                      ### find another peak which has the same value of the peak nearest to the center
		if abs(peaky[i]-peaky[dis_index])<tolerance:
			if abs(peakx[i]-peakx[dis_index])<min_dis and i!=dis_index and abs(y_intp[peakindex[i]-int(np/50)]-y_intp[peakindex[dis_index]-int(np/50)])<tolerance:
				min_dis=abs(peakx[i]-peakx[dis_index])
				min_index=i
  	if min_index==-9999999:
		print "\nWARNING!!!!!!!!!!!!!!!!!!!!!!!!!\nCan not find the range for averaging automatically AGAIN!!!\n PLEASE set the range manually!!!\nWARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		min_index=dis_index-1
	 	min_dis=1.0
		find_peak_flag=0
	else:
		print "\nWARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nSEEMS FIND THE PERIOD AFTER INCREASE THE TOLERANCE, CHECK THE PERIOD BY YOUR SELF!!!!\nWARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"

  num_ave=int(abs(peakindex[min_index]-peakindex[dis_index]))      ### determine the # of points for averaging
  if num_set!=0:
	num_ave=num_set
#  print peakx
#  print peakindex
  print "num_ave",num_ave
  if num_ave%2==0:
	spanl=num_ave/2-1
	spanr=num_ave/2
  else:
	spanl=num_ave/2
	spanr=num_ave/2
  x_ave=[]
  y_ave=[]
  for i in range(0,len(xd_intp)-spanr):                          ### plot the averaged curve
	sumy=0.0
	for j in range(i-spanl,i+spanr):
		sumy=sumy+y_intp[j]
	x_ave.append(xd_intp[i])
	y_ave.append(sumy/float(num_ave))
		
  dis_center_vac=100000.0 
  dis_vac_index=-999999
  for i in range(len(xd_intp)-spanr):                                  ### find the vacuum center
	dis_vac=abs(xd_intp[i]-vac_center)
	if dis_vac<dis_center_vac:
		dis_center_vac=dis_vac
		dis_vac_index=i
#  print "DIS_VAC"
#  print dis_center_vac,dis_vac_index,peakx[dis_vac_index]
#  print y_ave
 # plt.plot(xd_moved, vd0, 'o', xd_intp, y_intp, '-')
 # plt.plot(peakx, peaky, 'o', xd_intp, y_intp, '-')
#  xd,vd,vd0_intp = f_Smooth_Spline(xd0,vd0,width,ld,np)

#  print "V(z=0  )= %12.4f "%y_intp[0]  
#  print "V(z=1/2)= %12.4f "%y_intp[np/2] 

  print "\n ========================== RESULTS ===========================             \n"
  print "All the peaks and their index (check the average range is correct!!!)"
  peakx_novac=[]
  peakindex_novac=[]
#  for i in range(len(peakx)):
#	if peakx[i]<=zmax+0.5 and peakx[i]>=zmin-0.5:
#  		peakx_novac.append(peakx)	
  for i in range(len(peakindex)):
	if xd_intp[peakindex[i]]<=zmax+0.5 and xd_intp[peakindex[i]]>=zmin-0.5:
  		peakindex_novac.append(peakindex[i])
		peakx_novac.append(xd_intp[peakindex[i]])	
#  print peakx
#  print peakindex

#  print peakx_novac
  print '%16s %16s %16s' %('# of peak','correspond index','x coordinate')
  for i in range(len(peakx_novac)):
	print '%16d %16d %16.4f' %(i,peakindex_novac[i],peakx_novac[i])
  print "\n"
  print "current average peak and index:"
  print  peakindex[min_index]," ",xd_intp[peakindex[min_index]]
  print  peakindex[dis_index]," ",xd_intp[peakindex[dis_index]]
#  print peakindex_novac
  print "\n"
  print "1. Averaged Electric Potential at the Center of the Slab:"
  print y_ave[cent_index] 
  print "\n"
  print "2. The Coordinate of the Center of the Slab:"
  print x_ave[cent_index]
  print "\n"
  print "3. The Level of the Vacuum:"
  print y_ave[dis_vac_index]
  print "\n"
  print "4. E_vac-E_ave:"
  print y_ave[dis_vac_index]-y_ave[cent_index]
  print "\n"
  #print peakindex[dis_index],peakindex[min_index]
  
  print "Smoothed data and Averaged data is written to ",fname+".dat"
  ofile = open( fname+".dat",'w') 
  for i in range(len(xd_intp)):
    ofile.write('%12.6f %12.6f \n' %(xd_intp[i],y_intp[i]))
  ofile.write("\n")
  for i in range(len(x_ave)):
    ofile.write('%12.6f %12.6f \n' %(x_ave[i],y_ave[i]))
  ofile.close()

  savepath=os.getcwd()
  if find_peak_flag==1 or num_set!=0:
	plt.plot(peakx[dis_index], peaky[dis_index], 'o',peakx[min_index], peaky[min_index], 'o', xd_intp, y_intp, '-',x_ave,y_ave,'-',x_ave[cent_index],y_ave[cent_index],'o',x_ave[dis_vac_index],y_ave[dis_vac_index],'o')
	if save==1:
		plt.savefig(fname+'average.png')  
	else:
		plt.show()
  else:
	plt.plot(peakx, peaky, 'o',xd_intp, y_intp, '-',x_ave,y_ave,'-')
	if save==1:
		plt.savefig(fname+'average.png')  
	else:
		plt.show()
