#!/usr/env/bin python
import os
import numpy as np
import platform
from numpy.ctypeslib import ndpointer
import ctypes
from ctypes import byref, cdll, c_int, c_double,POINTER,CDLL,Structure,pointer

class c_complex128(Structure):
    _fields_ = [("real",c_double),("imag",c_double)]

    def to_complex(self):
        return self.real + 1j* self.imag

pt = ndpointer(dtype=np.complex128)
#pt = POINTER(c_double)

osname = platform.system()
if ("Windows" in osname):
#Default result with gfortran.exe -shared
    libpade = cdll.LoadLibrary(os.path.join(os.path.dirname(__file__),"./libpade.dll"))
    f_pade_multi = libpade.__m_pade_MOD_pade_multi 
    f_pade_test = libpade.__m_pade_MOD_pade_test 
    f_pade_test2 = libpade.__m_pade_MOD_pade_test_complex 
else:
    libpade = cdll.LoadLibrary(os.path.join(os.path.dirname(__file__),"./libpade.so.1")) 
    f_pade_multi = libpade.m_pade_mp_pade_multi_
    f_pade_test = libpade.m_pade_mp_pade_test_
    f_pade_test2 = libpade.m_pade_mp_pade_test_complex_
f_pade_multi.restype= None
#f_pade_multi.argtypes = [c_int,POINTER(c_double),POINTER(c_double),c_int,POINTER(c_double)
f_pade_multi.argtypes = [POINTER(c_int),pt,pt,POINTER(c_int),pt,pt,pt]


f_pade_test.argtypes = [POINTER(c_int),ndpointer(dtype=np.double)]
f_pade_test.restype = None


f_pade_test2.argtypes = [POINTER(c_int),POINTER(c_double)]
f_pade_test2.restype = None


#Load ac library
#This has dependancy
libac = CDLL(os.path.join(os.path.dirname(__file__),"./libacfreq.so"),mode=ctypes.RTLD_GLOBAL)
f_ac = libac.calcacfreq_
f_ac.argtypes = [POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
        ndpointer(dtype=np.double),
        ndpointer(dtype=np.complex128),
        POINTER(c_int),
        ndpointer(dtype=np.complex128),
        POINTER(c_complex128),
        POINTER(c_complex128),
        POINTER(c_complex128),
        ]
f_ac.restype = None


def get_array_p(x,typev):
    '''
    Convert a list to ndarray with specific type
    '''
    if (isinstance(x,list) or isinstance(x,tuple) or  x.dtype != typev):
        x2 = np.array(x,dtype=typev)
        return x2
    else:
        return x


def get_dcomplex_p(x):
    '''
    Convert a list to complex ndarray 
    Return ndarray
    :param x: list or ndarray
    '''
    if (isinstance(x,list) or isinstance(x,tuple) or  x.dtype != np.complex128):
        x2 = np.array(x,dtype=np.complex128)
        return x2
    else:
        return x
#   return x.ctypes.data_as(POINTER(c_double))

def pade_test(z):
    '''
    Test function
    '''
    n = len(z)
    n1 = c_int(n)
    print(n)
    f_pade_test(byref(n1),z)
    return z

def pade_test2(z):
    '''
    Test function
    '''
    n = len(z)
    n1 = c_int(n)
#   f_pade_test2(byref(n1),get_dcomplex_p(z))
    f_pade_test2(byref(n1),z.ctypes.data_as(POINTER(c_double)))
    return z
    

def pade_multi(z,f,zz):
    '''
    :param z: input of the function
    :param f: output of the function
    :param zz: where we need to Pade approximate to
    :return: two list, one with values and the other with derivatives
    '''
    n = len(z) #Complex 
    n2 = len(zz)
    z1 = get_dcomplex_p(z)
    f1 = get_dcomplex_p(f)
    zz1 = get_dcomplex_p(zz)
    pade_v = np.empty(n2,dtype=np.complex128)
    pade_d =  np.empty(n2,dtype=np.complex128)
    pade_v1= get_dcomplex_p(pade_v)
    pade_d1 = get_dcomplex_p(pade_d)
    n1 = c_int(n)
    n21 = c_int(n2)
#   print(f_pade_multi.argtypes[1])
    f_pade_multi(byref(n1),z1,f1,byref(n21),zz1,pade_v1,pade_d1)
#   f_pade_multi(n,z,f,n2,zz,pade_v,pade_d)
    return pade_v,pade_d

def ac_single(iac,omega,sc,npar,z):
    '''
    Analytical continuation, by Pade or multipole
    to a single complex position
    Refer to gap2a src_acfreq to see the meaning
    '''
    init = 0
    z1 = c_complex128(z,0)
    fz = c_complex128(1.0,1.0)
    dfz = c_complex128(1.0,1.0)

    nomega = len(omega)
    if (len(sc) != nomega):
        raise IndexError("Inconsistent array: omega(%i) and sc(%i)" % (nomega,len(sc)))

    if (nomega < npar):
        raise IndexError("npar=%i must be not larger than nomega=%i"%(npar,nomega))

    apar = np.empty(npar,dtype=np.complex128)

    omega1 = get_array_p(omega,np.double)
    sc1 = get_array_p(sc,np.complex128)
    apar1 = get_array_p(apar,np.complex128)

    f_ac(byref(c_int(init)),byref(c_int(iac)),byref(c_int(nomega)),
        omega1,
        sc1,
        byref(c_int(npar)),
        apar,
        pointer(z1),
        pointer(fz),pointer(dfz))

    return fz.to_complex(),dfz.to_complex()
            

def ac_multi(iac,npar,omega,sc,z):
    '''
    invoke ac_single to multiple values

    :param iac: method, 0=pade, 1=multi poles
    :param omega: positions on the imaginary axis to fit, must be real values (used as omega*1i)
    :param sc: a list, values to fit
    :param npar: the number of parameters in pade or mp, must be less than the number of input positions
    :param z: lists of complex where AC will be done
    :return: two lists, one values and one derivatives
    '''
    result =  zip(*[ac_single(iac,omega,sc,npar,x) for x in z])
    if (len(result) == 0):
        return [],[]
    else:
        return result

def ac_split(freqim,v,freqre,funcac):
    '''
    AC imag+ to real+ and imag- to real-, when 0 is a pole
    
    :param freqim: positions on the positive imaginary axis
    :param v: values
    :param freqre: position on the real axis where new values will be fit
    :param funcac: a function do AC, accept x/y/z and return y2,dy2 (ndarray or list)
    '''
#   freqimpos = [1j*x[0] for x in datai]
#   vimpos = [x[1]+1j*x[2] for x in datai]
#   freqimneg = [-1j*x[0] for x in datai]
#   vimneg = [x[1]-1j*x[2] for x in datai]
    freqimpos = freqim
    freqimneg = [-x for x in freqim]
    vimpos = v
    vimneg = [x.conjugate() for x in vimpos]
    result = [  funcac(freqimpos,vimpos,[x]) if x > 0  else funcac(freqimneg,vimneg,[x]) for x in freqre  ]
#Take away the first level array
    result = [ [x[0][0],x[1][0]] for x in result]
    return zip(*result)


