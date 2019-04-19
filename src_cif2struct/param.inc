
      parameter (mxsym=192,mxcen=8)
      common/comsgt/  iga(454),ipga(454),idla(454),itxta(454), &
                      grupaa(454)
      common/cellc/   Cell(6)
      common/symm/  ns,ncs,NcSym,nvt,mr6(36,mxsym),rm6(36,mxsym), &
                    s6(6,mxsym),vt6(6,mxcen),centr,latt,ngrupa, &
                    CrSystem,Monoclinic,BratSym(mxsym), &
                    ndim,ndimi,ndimq,ndimq1,ShiftSg(6), &
                    lattice,grupa,symmc(6,mxsym),smbx,smbc
      common/units/ dta,out,FormA1,FormA80,Cislo
      character*30 itxta,grupa
      character*15 cislo,symmc
      character*8  grupaa,smbc
      character*7  FormA1
      character*5  FormA80
      character*4  lattice
      character*1  smbx(3)
      integer dta,out,CrSystem
      logical NcSym
