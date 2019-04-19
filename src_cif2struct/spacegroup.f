      subroutine spacegroup

      use struct, only : sgname,aa,bb,cc,alpha,&
           pos,mult,aname,nat

      implicit real*8 (a-h,o-z)
      include 'param.inc'
      dimension x(3,999),y(3),xx(3),yy(3),zz(3),zzr(3),g(3,3)
      character*80 Radka
      character*8  Atom
      character*5  t5
      character*2  nty
      real origin(3)

      data torad/0.017453293/

      read(1,*) (origin(i),i=1,3)

      indatom=0 
1000  read(1,*) (cell(i),i=1,6)
!      write(*,'(6f10.6)') cell
      do 1100 i=1,3
        g(i,i)=cell(i)**2
1100  continue
      g(1,2)=Cell(1)*Cell(2)*cos(torad*Cell(6))
      g(1,3)=Cell(1)*Cell(3)*cos(torad*Cell(5))
      g(2,3)=Cell(2)*Cell(3)*cos(torad*Cell(4))
      g(2,1)=g(1,2)
      g(3,1)=g(1,3)
      g(3,2)=g(2,3)
2000  read(1,*) Grupa
!      write(*,'(a)') grupa
      call ChangeToTitle(Grupa)
      call PitSym(ich)
      if(ich.ne.0.) then
         write(*,'(2a)') 'wrong Grupa: ',Grupa
         stop         
      endif  
      n=0
      do 2400 i=1,ns
        do 2350 k=1,ndim
          n=max0(n,idel(symmc(k,i))+1)
2350    continue
2400  continue
      if(n.lt.5) n=5
      do 2500 i=1,ns
        radka=' '
        k=1
        do 2450 j=1,ndim
          radka(k+n-idel(symmc(j,i)):)=symmc(j,i)
          k=k+n
2450    continue
2500  continue
      na=0
3000  na=na+1
      call zhusti(t5)
      read(1,*,end=9000) atom
      aname(na)=atom
!      write(*,'(a)') atom
      if(Atom.eq.' ') go to 9000
      read(1,*) (xx(i),i=1,3) 
      xx(1:3)=xx(1:3)+origin(1:3)
!      write(*,'(3f10.7)') xx
      n=0
      do 4000 i=1,ns
        call multm(rm6(1,i),xx,y,3,3,1)
        do 3500 j=1,3-ncs
          if(j.eq.2) then
            do 3120 m=1,3
              y(m)=-y(m)
3120        continue
          endif
          do 3110 m=1,3
            yy(m)=y(m)+s6(m,i)
3110      continue
          call od0do1(yy,yy,3)
          do 3300 k=1,n
            do 3200 ivt=1,nvt
              do 3150 m=1,3
                zz(m)=yy(m)+vt6(m,ivt)-x(m,k)
                zz(m)=zz(m)-anint(zz(m))
3150          continue
              call multm(g,zz,zzr,3,3,1)
              dist=sqrt(scalmul(zz,zzr))
              if(dist.lt..00001d0) go to 3500
3200        continue
3300      continue
          n=n+1
          call CopyVek(yy,x(1,n),3)
3500    continue
4000  continue
      Radka=Atom
      k=idel(Radka)+1
      do 5000 i=1,n
        write(Radka(k:),'(''/'',i3)') i
        call Zhusti(Radka)
        call ChangeToTitle(Radka)
!        write(out,'(1x,a10,3f12.8)') Radka,(x(j,i),j=1,3)
        indatom=indatom+1
        pos(1:3,indatom)=x(1:3,i)         
5000  continue
      mult(na)=n 
      go to 3000
9000  continue

      sgname=Grupa       
      aa=cell(1)
      bb=cell(2)
      cc=cell(3)
      alpha(1)=cell(4)
      alpha(2)=cell(5)
      alpha(3)=cell(6)
      nat=na-1

      end

      function idel(l)
      character*(*) l
      idel=len(l)
1000  if(idel.ne.0) then
        if(l(idel:idel).eq.' ') then
          idel=idel-1
          go to 1000
        endif
      endif
      return
      end
      subroutine ChangeToTitle(at)
      implicit real*8 (a-h,o-z)
      character*(*) at
      call mala(at)
      call velka(at(1:1))
      return
      end


      subroutine velka(veta)
      implicit real*8 (a-h,o-z)
      character*(*) veta
      character*1 malchar(26),velchar(26),v
      logical first,souvisle
      save na,nz,ndif
      data first,souvisle/2*.true./
      data malchar/'a','b','c','d','e','f','g','h','i','j','k','l','m', &
                   'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      data velchar/'A','B','C','D','E','F','G','H','I','J','K','L','M', &
                   'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      if(first) then
        do 1000 i=1,25
          if(ichar(malchar(i+1))-ichar(malchar(i)).ne.1) go to 1100
          if(ichar(velchar(i+1))-ichar(velchar(i)).ne.1) go to 1100
1000    continue
        na=ichar('a')
        nz=ichar('z')
        ndif=ichar('A')-na
        go to 1200
1100    souvisle=.false.
1200    first=.false.
      endif
      if(souvisle) then
        do 2000 i=1,idel(veta)
         n=ichar(veta(i:i))
         if(n.ge.na.and.n.le.nz) veta(i:i)=char(n+ndif)
2000    continue
      else
        do 4000 i=1,idel(veta)
          v=veta(i:i)
          do 3000 j=1,26
            if(v.eq.malchar(j)) go to 3500
3000      continue
          go to 4000
3500      veta(i:i)=velchar(j)
4000    continue
      endif
      return
      end
      subroutine mala(veta)
      implicit real*8 (a-h,o-z)
      character*(*) veta
      character*1 malchar(26),velchar(26),v
      logical first,souvisle
      save na,nz,ndif
      data first,souvisle/2*.true./
      data malchar/'a','b','c','d','e','f','g','h','i','j','k','l','m', &
                   'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      data velchar/'A','B','C','D','E','F','G','H','I','J','K','L','M', &
                   'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      if(first) then
        do 1000 i=1,25
          if(ichar(malchar(i+1))-ichar(malchar(i)).ne.1) go to 1100
          if(ichar(velchar(i+1))-ichar(velchar(i)).ne.1) go to 1100
1000    continue
        na=ichar('A')
        nz=ichar('Z')
        ndif=ichar('a')-na
        go to 1200
1100    souvisle=.false.
1200    first=.false.
      endif
      if(souvisle) then
        do 2000 i=1,idel(veta)
        n=ichar(veta(i:i))
        if(n.ge.na.and.n.le.nz) veta(i:i)=char(n+ndif)
2000    continue
      else
        do 4000 i=1,idel(veta)
          v=veta(i:i)
          do 3000 j=1,26
            if(v.eq.velchar(j)) go to 3500
3000      continue
          go to 4000
3500      veta(i:i)=malchar(j)
4000    continue
      endif
      return
      end
      subroutine PitSym(ich)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      dimension tau(4)
      character*80 t80,Veta
      character*30 itxt
      logical BylaNeni
      MonoclinicDef=Monoclinic
      if(MonoclinicDef.lt.1.or.MonoclinicDef.gt.3) MonoclinicDef=2
      call CheckSystem(Cell,Monoclinic,CrSystem)
      call Zhusti(Grupa)
1100  if(Grupa(1:1).eq.'-') then
        call ChangeToTitle(Grupa(2:))
      else
        call ChangeToTitle(Grupa)
      endif
      i=index(Grupa,';')
      j=index(Grupa(2:),'x')
      k=index(Grupa,'y')
      l=index(Grupa,'z')
      if(i.le.0.and.j.le.0.and.k.le.0.and.l.le.0.and.grupa(1:1).ne.'-') &
         then
        Lattice=Grupa(1:1)
        if(Lattice.eq.'X') Grupa(1:1)='P'
        imd=0
        BylaNeni=.false.
        do 2000 i=1,454
          if(imd.ne.0) imd=imd+1
          if((Grupaa(i).eq.'P2'.and.imd.eq.0).or.Grupaa(i).eq.'Pmm2') &
             imd=1
          if(Grupa.eq.Grupaa(i)) then
            if(ipga(i).ge.3.and.ipga(i).le.5) then
              if(mod(imd-1,3)+1.ne.Monoclinic) then
                BylaNeni=.true.
                go to 2000
              endif
            endif
            go to 2100
          endif
2000    continue
        if(BylaNeni) then
          go to 8000
        else
          go to 8170
        endif
2100    k=index(itxta(i),'#')
        if(k.gt.0) then
          itxt=itxta(i)(1:k-1)
          Veta=itxta(i)(k+1:idel(itxta(i)))
          do 2105 j=1,idel(Veta)
            if(Veta(j:j).eq.',') Veta(j:j)=' '
2105      continue
          call StToReal(Veta,ShiftSg,3,.false.,ich)
        else
          itxt=itxta(i)
          do 2120 j=1,ndim
            ShiftSg(j)=0.
2120      continue
        endif
        ipg=ipga(i)
        ngrupa=iga(i)
        idl=idla(i)
        if(ipg.le.2) then
          CrSystem=1
        else if(ipg.ge.3.and.ipg.le.5) then
          imd=mod(imd-1,3)+1
          if(CrSystem.eq.1) go to 8000
          CrSystem=2+Monoclinic*10
        else if(ipg.ge.6.and.ipg.le.8) then
          imd=mod(imd+1,3)+1
          if(mod(CrSystem,10).le.2.or.CrSystem.eq.6) go to 8000
          CrSystem=3
        else if(ipg.ge.9.and.ipg.le.15) then
          if(mod(CrSystem,10).le.3.or.CrSystem.eq.6) go to 8000
          CrSystem=4
        else if(ipg.ge.16.and.ipg.le.27) then
          if(CrSystem.ne.6) go to 8000
          if(ipg.le.20) CrSystem=5
        else if(ipg.ge.28.and.ipg.le.32) then
          if(CrSystem.ne.7) go to 8000
        endif
        if(mod(CrSystem,10).ne.2) Monoclinic=0
        Grupa(1:1)=Lattice
        if(itxt(1:1).eq.'-') then
          itxt(2:2)=Lattice
        else
          itxt(1:1)=Lattice
        endif
      else
        itxt=Grupa
        idl=1
        do 2150 i=1,idel(Grupa)
          if(Grupa(i:i).eq.';') idl=idl+1
2150    continue
        do 2160 j=1,ndim
          ShiftSg(j)=0.
2160    continue
      endif
      call GenCentr(itxt,nc,i1xx)
       if(i1xx.eq.1) goto 8110
       if(i1xx.eq.2) goto 9900
      call GenSym(itxt,tau,idl,*9900)
!5000  if(StdSg.ne.1) call OriginShift(ShSg)
5000  do 5400 i=1,ns
        calL CodeSym(rm6(1,i),s6(1,i),symmc(1,i))
        do 5300 j=1,ndim
          if(index(symmc(j,i),'?').gt.0) then
            ns=0
            go to 8140
          endif
5300    continue
5400  continue
      ich=0
      go to 9999
8000  call Chybne('cell parameters are not consistent with '// &
                  'space group','try another space group or '// &
                  'change cell parameters',0)
      go to 9900
8110  t80='the cell centring symbol is not correct'
      go to 8500
8140  t80='origin shift is not acceptable'
      go to 8500
8170  t80='the symbol wasn''t found on the list'
8500  Veta='incorrect space group symbol'
      call Chybne(Veta,t80,0)
9900  ich=1
9999  return
      end
      subroutine zhusti(s)
      implicit real*8 (a-h,o-z)
      character*(*) s
      idl=idel(s)
      k=0
      do 1000 i=1,idl
        if(s(i:i).eq.' ') go to 1000
        k=k+1
        if(k.ne.i) then
          s(k:k)=s(i:i)
          s(i:i)=' '
        endif
1000  continue
      return
      end
      subroutine StToReal(Veta,X,n,Message,ErrFlg)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      dimension X(n)
      character*(*) Veta
      character*5 format
      logical Message
      integer ErrFlg
      ErrFlg=0
      if(idel(Veta).le.0) go to 9999
      k=0
      do 1000 i=1,n
        if(k.ge.len(Veta)) go to 9100
        call Kus(Veta,k,Cislo)
        j=index(Cislo,'/')
        if(j.le.0) then
          call posun(Cislo,1)
          read(Cislo,'(f15.0)',err=9000) X(i)
        else
          idl=idel(Cislo)
          format='(i  )'
          if(j.eq.1.or.j.eq.idl) then
            X(i)=0.
            go to 1000
          endif
          write(format(3:4),100) j-1
          read(Cislo(1:j-1),format,err=9000) i1
          write(format(3:4),100) idl-j
          read(Cislo(j+1:idl),format,err=9000) i2
          X(i)=float(i1)/float(i2)
        endif
1000  continue
      go to 9999
9000  errflg=1
      if(Message) then
        call Chybne('illegal real number - try again',Cislo,0)
      endif
      go to 9999
9100  errflg=1
      if(Message) then
        call Chybne('the string doesn''t contains all '// &
                    'requested real numbers',Veta(1:idel(Veta)),0)
      endif
9999  return
100   format(i2)
      end
      subroutine chybne(text1,text2,konec)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      character*132 ven
      character*(*) text1,text2
      if(konec.eq.1) then
        ven=' Fatal -'
      else if(konec.eq.0) then
        ven=' Error -'
      else
        ven=' Warning -'
      endif
      ven=ven(1:idel(ven))//' '//text1(1:idel(text1))
      n=idel(ven)
      write(out,FormA1)(ven(i:i),i=1,min0(79,n))
      n=idel(text2)
      if(n.ne.0) then
        if(konec.eq.-1) then
          ven='           '//text2(1:idel(text2))
        else
          ven='         '//text2(1:idel(text2))
        endif
        n=idel(ven)
        write(out,FormA1)(ven(i:i),i=1,min0(79,n))
      endif
      if(konec.eq.1) stop
9999  return
      end
      subroutine kus(radka,k,slovo)
      implicit real*8 (a-h,o-z)
      character*(*) radka,slovo
      mxs=len(slovo)
      mxr=len(radka)
      slovo=' '
      l=0
      if(k.ge.mxr) go to 9999
      k=k+1
1000  if(radka(k:k).eq.' ') then
        if(k.ge.mxr) go to 9999
        k=k+1
        go to 1000
      endif
2000  if(radka(k:k).ne.' ') then
        if(l.lt.mxs) then
          l=l+1
          slovo(l:l)=radka(k:k)
        endif
        if(k.ge.mxr) go to 9999
        k=k+1
        go to 2000
      endif
3000  if(radka(k:k).eq.' ') then
        if(k.ge.mxr) go to 9999
        k=k+1
        go to 3000
      endif
      k=k-1
9999  return
      end
      subroutine posun(cislo,tecka)
      implicit real*8 (a-h,o-z)
      character*(*) cislo
      integer tecka
      l=idel(cislo)
      m=len(cislo)
      if(tecka.eq.1.and.index(cislo,'.').le.0.and.l.lt.m) then
        l=l+1
        cislo(l:l)='.'
      endif
      if(l.lt.m) then
        do 1000 i=l,1,-1
          cislo(m-l+i:m-l+i)=cislo(i:i)
          cislo(i:i)=' '
1000    continue
      endif
      return
      end
      subroutine CheckSystem(Cell,Monoclinic,CrSystem)
      implicit real*8 (a-h,o-z)
      dimension Cell(6)
      integer CrSystem
      j=0
      Monoclinic=0
      do 1000 i=1,3
        if(abs(Cell(i+3)-90.).lt..001d0) then
          j=j+1
        else
          Monoclinic=i
        endif
1000  continue
      if(j.ne.2) then
        Monoclinic=0
      else
        if(Monoclinic.ne.0) then
          if(abs(Cell(Monoclinic+3)-120.).lt..001d0) then
            CrSystem=6
          else
            CrSystem=Monoclinic*10+2
          endif
        endif
        go to 2000
      endif
      CrSystem=1
      if(j.eq.3) then
        if(abs(Cell(1)-Cell(2)).gt..0001d0) then
          CrSystem=3
        else if(abs(Cell(1)-Cell(3)).gt..0001d0) then
          CrSystem=4
        else
          CrSystem=7
        endif
      else if(j.eq.0) then
        do 1100 i=2,3
          if(abs(Cell(1)-Cell(i)).gt..0001d0) go to 2000
          if(abs(Cell(4)-Cell(i+3)).gt..0001d0) go to 2000
1100    continue
        CrSystem=5
      endif
2000  return
      end
      subroutine CodeSym(rmm,s,is)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      character*15 is(6)
      dimension rmm(36),s(6)
      do 5000 i=1,ndim
        ij=i
        is(i)=' '
        do 1000 k=1,100
          p=s(i)*float(k)
          l=nint(p)
          if(abs(p-float(l)).lt..0001d0) go to 1100
1000    continue
        is(i)='?/?'
        go to 1200
1100    if(l.ne.0) then
          write(is(i)(1:8),'(i4,''/'',i3)') l,k
          if(k.eq.1) is(i)(5:8)='  '
        endif
1200    call zhusti(is(i))
        kk=idel(is(i))
        do 4000 j=1,ndim
          l=nint(rmm(ij))
          if(iabs(l).gt.0) then
            if(l.gt.0.and.kk.gt.0) then
              kk=kk+1
              is(i)(kk:kk)='+'
            else if(l.lt.0) then
              kk=kk+1
              is(i)(kk:kk)='-'
            endif
            if(iabs(l).gt.1) then
              write(is(i)(kk+1:),'(i2)') iabs(l)
            endif
            call zhusti(is(i))
            kk=max0(idel(is(i)),1)
            is(i)=is(i)(1:kk)//smbx(j)
            kk=kk+1
            call zhusti(is(i))
            kk=idel(is(i))
          endif
          ij=ij+ndim
4000    continue
5000  continue
      return
      end
      subroutine GenCentr(itxt,nc,i1xx)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      character*(*) itxt
      if(itxt(1:1).eq.'-') then
        ncs=1
        itxt=itxt(2:)
      else if(itxt(1:1).eq.'+') then
        ncs=2
        itxt=itxt(2:)
      else
        ncs=2
      endif
      NcSym=ncs.eq.2
      nc=index(smbc,itxt(1:1))
      if(nc.gt.0.and.nc.le.8) then
        call GenVecCentr(nc,ich)
        if(ich.ne.0) then
        i1xx=2
        return
        endif
      else
        i1xx=1
        return
      endif
      Lattice=itxt(1:1)
      itxt=itxt(2:)
      i1xx=0
      return
      end
      subroutine GenVecCentr(nc,ich)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      ich=0
      if(nc.lt.8) then
        nvt=1
        do 1005 i=1,ndim
          do 1000 j=1,4
            vt6(i,j)=0.
1000      continue
1005    continue
      endif
      if(nc.gt.1.and.nc.le.5) then
        do 1021 i=1,3
          vt6(i,2)=.5
1021    continue
        if(nc.ne.5) vt6(nc-1,2)=0.
        nvt=2
      else if(nc.eq.6) then
        nvt=3
        vt6(1,2)=.3333333333333d0
        vt6(2,2)=.6666666666667d0
        vt6(3,2)=.6666666666667d0
        vt6(1,3)=.6666666666667d0
        vt6(2,3)=.3333333333333d0
        vt6(3,3)=.3333333333333d0
      else if(nc.eq.7) then
        nvt=4
        do 1042 j=2,4
          do 1041 i=1,3
            if(i.ne.j-1) then
              vt6(i,j)=.5
            else
              vt6(i,j)=0.
            endif
1041      continue
1042    continue
      endif
      return
      end
      subroutine GenSym(itxt,tau,idl,*)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      dimension tau(idl),vtp(6),tv(3,8),tvo(3,6),pp(6)
      character*(*) itxt
      character*2 nty
      character*5 ir
      character*20 ito
      character*80 t80,errtxt
      data ir/'12346'/,ito/'xyz"''*12345abcnuvwd;'/
      data tv/.5,.0,.0,.0,.5,.0,.0,.0,.5,.5,.5,.5,.25,.0,.0, &
              .0,.25,.0,.0,.0,.25,.25,.25,.25/
      data tvo/1.,0.,0.,0.,1.,0.,0.,0.,1.,1.,1.,0.,1.,-1.,0.,1.,1.,1./
      ns=1
      call UnitMat(rm6(1,1),ndim)
      k=0
      do 1100 i=1,ndim
        do 1000 j=1,ndim
          k=k+1
          mr6(k,1)=nint(rm6(k,1))
1000    continue
        vt6(i,1)=0.
1100  continue
      id=idel(itxt)
2000  if(id.le.0) go to 9999
      if(itxt(1:1).eq.'-') then
        iz=-1
        itxt=itxt(2:)
        id=id-1
      else if(itxt(1:1).eq.'+') then
        iz= 1
        itxt=itxt(2:)
        id=id-1
      else
        iz= 1
      endif
      if(id.le.0) then
        errtxt='rotational information is missing'
        go to 8000
      endif
      nr=index(ir,itxt(1:1))
      itxt=itxt(2:)
      if(nr.le.0) then
        errtxt='incorrect rotational part'
        go to 8000
      else if(nr.eq.1) then
        if(idel(itxt).gt.0) then
          errtxt='incorrect rotational part'
          go to 8000
        else
          go to 5000
        endif
      endif
      ior=-1
      do 2100 i=1,ndim
        vtp(i)=0.
2100  continue
      pom=0.
3000  if(idel(itxt).gt.0) then
        i=index(ito,itxt(1:1))
        if(i.gt.0) then
          itxt=itxt(2:)
        else
          errtxt='incorrect orientation or translation part'
          go to 8000
        endif
      else
        i=20
      endif
      if(i.le.6) then
        if(ior.eq.-1) then
          ior=i
          go to 3000
        else
          errtxt='orientation is doubled'
          go to 8000
        endif
      else if(i.le.19) then
        nt=i-6
        if(nt.le.5) then
          p=nr
          if(nr.eq.5) then
            p=6.
          else if(nr.eq.6.or.nr.eq.7) then
            p=2.
          else if(nr.eq.8) then
            p=3.
          endif
          pom=pom+float(nt)/p
        else
          do 3100 i=1,3
            vtp(i)=vtp(i)+tv(i,nt-5)
3100      continue
        endif
      else if(i.eq.20) then
        if(ior.lt.0) then
          if(ns.eq.1) then
            ior=3
          else if(ns.eq.2.or.ns.eq.3) then
            if(nr.eq.2.and.(nrold.eq.2.or.nrold.eq.4)) ior=1
            if(nr.eq.2.and.(nrold.eq.3.or.nrold.eq.5)) ior=5
            if(nr.eq.3) ior=6
          endif
        endif
        if(ior.lt.0) then
          errtxt='orientation is not defined'
          go to 8000
        endif
        if(ior.eq.6) then
          if(nr.eq.3) then
            nr=8
            ior=1
          else
            errtxt='operator cannot have this orientation'
            go to 8000
          endif
        else if(ior.eq.4.or.ior.eq.5) then
          if(nr.eq.2) then
            nr=ior+2
            if(ns.ne.1) then
              ior=iorold
            else
              ior=3
            endif
          else
            errtxt='operator cannot have this orientation'
            go to 8000
          endif
        endif
        ns=ns+1
        call SetGen(mr6(1,ns),nr-1,ior,ndim)
        call MatVek(mr6(1,ns),ShiftSg,pp,ndim)
        do 3200 i=1,ndim
          if(iz.eq.-1) then
            k=(i-1)*ndim
            do 3150 j=1,ndim
              ind=j+k
              mr6(ind,ns)=-mr6(ind,ns)
3150        continue
            pp(i)=-pp(i)
          endif
          if(i.le.3) then
            s6(i,ns)=vtp(i)+ShiftSg(i)-pp(i)
            s6(i,ns)=s6(i,ns)+pom*tvo(i,ior)
          endif
3200    continue
        call od0do1(s6(1,ns),s6(1,ns),ndim)
        call NormCentr(s6(1,ns))
        iorold=ior
        nrold=nr
        if(idel(itxt).gt.0) then
          go to 2000
        else
          go to 5000
        endif
      endif
      go to 3000
5000  call CompleteSym(ich)
      if(ich.ne.0) go to 8100
      do 5100 k=1,ns
        l=0
        do 5020 i=1,ndim
          do 5010 j=1,ndim
            l=l+1
            rm6(l,k)=mr6(l,k)
5010      continue
5020    continue
5100  continue
      go to 9999
8000  write(t80,'(''in the '',i1,a2,'' Hall''''s generator'')') &
            ns,nty(ns)
      call Chybne(t80,errtxt(1:idel(errtxt)),0)
8100  return 1
9999  return
      end
      subroutine SetGen(mr,nr,md,ndim)
      implicit real*8 (a-h,o-z)
      dimension irm(9,7),mr(*)
      data irm/-1,0,0,0,-1,0,0,0,1,0,1,0,-1,-1,0,0,0,1,0,1,0,-1,0,0, &
               0,0,1,1,1,0,-1,0,0,0,0,1,0,1,0,1,0,0,0,0,-1,0,-1,0,-1,0,0 &
              ,0,0,-1,0,1,0,0,0,1,1,0,0/
      id=md*3
      i1=-3
      mdd=md-1
      do 1000 i=1,ndim*ndim
        mr(i)=0
1000  continue
      do 2000 i=1,3
        i1=i1+3
        i1p=mod(i1+id,9)+1
        do 1500 j=1,3
          ij=i1+j
          ijp=i1p+mod(j+mdd,3)
          im=(ijp-1)/3+1
          jm=mod(ijp-1,3)+1
          mr(jm+(im-1)*ndim)=irm(ij,nr)
1500    continue
2000  continue
      return
      end
      subroutine CompleteSym(ich)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      dimension mrp(36),rmp(36),sp(6),mrpm(36),is(mxsym)
      logical Eqrv,Eqiv
      ich=0
      nss=ns
      do 400 i=1,mxsym
        is(i)=0
400   continue
500   i=0
1000  if(i.ge.ns) go to 9000
      i=i+1
      j=is(i)
2000  if(j.ge.ns) go to 1000
      j=j+1
      call multmi(mr6(1,i),mr6(1,j),mrp,ndim,ndim,ndim)
      do 2100 k=1,ndimq
        rmp(k)=mrp(k)
2100  continue
      call MatVek(mr6(1,i),s6(1,j),sp,ndim)
      do 2300 k=1,ndim
        sp(k)=sp(k)+s6(k,i)
2300  continue
      call od0do1(sp,sp,ndim)
      call NormCentr(sp)
      do 5000 k=1,ns
        do 4000 l=1,3-ncs
          if(l.eq.1) then
            call CopyMatI(mrp,mrpm,ndim)
          else
            do 3000 m=1,ndimq
              mrpm(m)=-mrp(m)
3000        continue
          endif
          if(Eqiv(mrpm,mr6(1,k),ndimq)) then
            if(Eqrv(sp,s6(1,k),ndim,.0001d0)) then
              go to 2000
            else
              call Chybne('symmetry operators are not '// &
                          'consistent','with translation symmetry',0)
              go to 8000
            endif
          endif
4000    continue
5000  continue
      ns=ns+1
      if(ns.gt.mxsym) then
        call Chybne('the number of symmetry operators'// &
                    ' exceeds the limit','check used symbol',0)
        go to 8000
      else
        call CopyMatI(mrp,mr6(1,ns),ndim)
        call CopyMat (rmp,rm6(1,ns),ndim)
        call CopyVek(sp,s6(1,ns),ndim)
        call CodeSym(rm6(1,ns),s6(1,ns),symmc(1,ns))
      endif
      is(i)=j
      go to 500
8000  ns=nss
      ich=1
9000  return
      end
      subroutine od0do1(x,y,n)
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n)
      do 1000 i=1,n
        y(i)=x(i)-aint(x(i))
        if(y(i).lt.0.) y(i)=y(i)+1.d0
        if(y(i).gt..999999.or.y(i).lt..000001d0) y(i)=0.d0
1000  continue
      return
      end
      subroutine MatVek(m,v,p,n)
      implicit real*8 (a-h,o-z)
      dimension m(n,n),v(n),p(n)
      do 2000 i=1,n
        pi=0.
        do 1000 j=1,n
          pi=pi+float(m(i,j))*v(j)
1000    continue
        p(i)=pi
2000  continue
      return
      end
      subroutine UnitMat(a,n)
      implicit real*8 (a-h,o-z)
      dimension a(n,n)
      do 2000 i=1,n
        do 1000 j=1,n
          if(i.eq.j) then
            a(i,j)=1.
          else
            a(i,j)=0.
          endif
1000    continue
2000  continue
      return
      end
      character*2 function nty(i)
      character*2 th(4)
      data th/'st','nd','rd','th'/
      j=min0(4,i)
      nty=th(j)
      return
      end
      subroutine NormCentr(x)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      dimension x(*)
      do 5000 i=1,ndim
        xp=100.
        n=0
        do 3000 j=1,nvt
          do 2500 k=1,i-1
            if(vt6(k,j).ne.0.) go to 3000
2500      continue
          if(vt6(i,j).eq.0.) go to 3000
          if(vt6(i,j).lt.xp) then
            n=j
            xp=vt6(i,j)
          endif
3000    continue
        if(n.eq.0) go to 5000
3100    if(x(i).ge.xp-.0001d0) then
          do 3150 k=i,ndim
            x(k)=x(k)-vt6(k,n)
            if(x(k).lt.-.0001d0.and.k.ne.i) x(k)=x(k)+1.d0
3150      continue
          go to 3100
        endif
3200    if(x(i).lt.-.0001d0) then
          do 3250 k=i,ndim
            x(k)=x(k)+vt6(k,n)
            if(x(k).ge.1.d0.and.k.ne.i) x(k)=x(k)-1.d0
3250      continue
          go to 3200
        endif
3300    if(abs(x(i)).lt..0001d0) x(i)=0.d0
5000  continue
      return
      end
      subroutine CopyMat(a,b,n)
      implicit real*8 (a-h,o-z)
      dimension a(n,n),b(n,n)
      do 2000 i=1,n
        do 1000 j=1,n
          b(i,j)=a(i,j)
1000    continue
2000  continue
      return
      end
      subroutine CopyMatI(a,b,n)
      implicit real*8 (a-h,o-z)
      integer a(n,n),b(n,n)
      do 2000 i=1,n
        do 1000 j=1,n
          b(i,j)=a(i,j)
1000    continue
2000  continue
      return
      end
      subroutine CopyVek(a,b,n)
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n)
      do 1000 i=1,n
        b(i)=a(i)
1000  continue
      return
      end
      logical function eqiv(a,b,n)
      integer a(n),b(n)
      eqiv=.false.
      do 1000 i=1,n
        if(a(i).ne.b(i)) go to 2000
1000  continue
      eqiv=.true.
2000  return
      end
      logical function eqrv(a,b,n,dif)
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n)
      eqrv=.false.
      do 1000 i=1,n
        if(abs(a(i)-b(i)).ge.dif) go to 9999
1000  continue
      eqrv=.true.
9999  return
      end
      subroutine multm(a,b,c,n1,n2,n3)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*),c(*)
      do 3000 i=1,n1
        jkp=1
        ik=i
        do 2000 k=1,n3
          ij=i
          jk=jkp
          p=0.
          do 1000 j=1,n2
            p=p+a(ij)*b(jk)
            ij=ij+n1
            jk=jk+1
1000      continue
          c(ik)=p
          ik=ik+n1
          jkp=jkp+n2
2000    continue
3000  continue
      return
      end
      subroutine multmi(a,b,c,n1,n2,n3)
      implicit real*8 (a-h,o-z)
      integer a(*),b(*),c(*)
      do 3000 i=1,n1
        jkp=1
        ik=i
        do 2000 k=1,n3
          ij=i
          jk=jkp
          ip=0
          do 1000 j=1,n2
            ip=ip+a(ij)*b(jk)
            ij=ij+n1
            jk=jk+1
1000      continue
          c(ik)=ip
          ik=ik+n1
          jkp=jkp+n2
2000    continue
3000  continue
      return
      end
     subroutine ReadText(text,stout,stin)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      character*(*) text,stout,stin
      character*120 format
      character*80 t80
      format=' : '
      i=idel(text)
      j=idel(stin)
      if(j.gt.0) format=' ['//stin(1:j)// ']'//format(1:3)
      if(i.gt.0) format=text(1:i)//format(1:idel(format)+1)
      if(i.le.0) format='?'//format(1:idel(format)+1)
      write(out,'('' '//format(1:idel(format)+1)//''',$)')
      read(dta,FormA80) t80
      if(idel(t80).le.0) then
        stout=stin
      else
        stout=t80
      endif
      return
      end
      subroutine ReadReal(text,a,da,n,iw,ierr)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      dimension a(n),da(n)
      character*(*) text
      character*120 String,Cisla
      character*80 veta
      character*8  fwd,fa1
      character*2  nty
      data fwd,fa1/'(f15.12)','(99a1,$)'/
      if(n.le.0) return
      np=1
      ierr=0
1000  String=' '//text(1:idel(text))
      call ShortApostr(String)
      ifk=idel(String)
      if(iw.gt.0) then
        String=String(1:ifk)//' ['
        Cisla=' '
        ip=1
        do 1050 i=np,n
          if(abs(a(i)).gt.0.) then
            j=max0(6-mantisa10(a(i)),0)
            j=min0(j,6)
          else
            j=6
          endif
          write(fwd(6:7),'(i2)') j
          write(Cisla(ip:),fwd) a(i)
          call ZdrcniCisla(Cisla,i)
          ip=idel(Cisla)+1
1050    continue
        String=String(1:idel(String))//Cisla(:ip-1)//'] : '
      else
        String=String(1:idel(String))//' : '
      endif
      ip=idel(String)+1
      write(fa1(2:3),'(i2)') ip
      write(out,fa1)(String(i:i),i=1,ip)
      read(dta,FormA80) veta
      idelv=idel(veta)
      if(idelv.le.0.and.iw.ne.0) then
        do 1100 i=np,n
          a(i)=da(i)
1100    continue
      else
        if(idelv.eq.0) then
          ierr=1
          return
        endif
        k=0
        do 1200 i=np,n
          call kus(veta,k,cislo)
          a(i)=ReadRealFract(cislo,ich)
          if(ich.ne.0) go to 1250
          if(k.eq.80) go to 1260
1200    continue
        return
1250    np=i
        write(out,'('' Error - non-numerical character in the '',i1,a2,''item of this line'')') np,nty(np)
        go to 1000
1260    np=i+1
        if(np.gt.n) return
        write(out,'('' Response not complete continue with the '',i1,a2,'' item of this line'')') np,nty(np)
        go to 1000
      endif
101   format(i2)
102   format(' String contains non-numerical character, try again')
      end
      function ReadRealFract(veta,ich)
      implicit real*8 (a-h,o-z)
      character*(*) veta
      character*15 pomc
      character*5 format
      ich=0
      pomc=veta
      idl=idel(pomc)
      i=index(pomc,'/')
      if(i.le.0) then
        call posun(pomc,1)
        read(pomc,'(f15.0)',err=9000) ReadRealFract
      else
        format='(i  )'
        if(i.eq.1.or.i.eq.idel(pomc)) then
          ReadRealFract=0.
          go to 9999
        endif
        write(format(3:4),100) i-1
        read(pomc(1:i-1),format,err=9000) i1
        write(format(3:4),100) idl-i
        read(pomc(i+1:idl),format,err=9000) i2
        ReadRealFract=float(i1)/float(i2)
      endif
      go to 9999
9000  ich=1
9999  return
100   format(i2)
      end
      subroutine ShortApostr(format)
      implicit real*8 (a-h,o-z)
      character*(*) format
      character*80  t80
      j=1
      t80=' '
      do 1000 i=1,idel(format)
        if(format(i:i).eq.'''') then
          t80(j:j+1)=''''''
          j=j+2
        else
          t80(j:j)=format(i:i)
          j=j+1
        endif
        if(j.eq.80) go to 2000
1000  continue
2000  format=t80(1:j)
      return
      end
      function Mantisa10(x)
      implicit real*8 (a-h,o-z)
      pom=abs(x)
      Mantisa10=0
1000  if(pom.gt.10.) then
        pom=pom*.1
        Mantisa10=Mantisa10+1
        go to 1000
      endif
2000  if(pom.le.1.) then
        pom=pom*10.
        Mantisa10=Mantisa10-1
        go to 2000
      endif
      return
      end
      subroutine ZdrcniCisla(radka,kolik)
      implicit real*8 (a-h,o-z)
      character*(*) radka
      character*80 t80
      character*15 t15
      if(idel(radka).eq.0) return
      k=0
      do 2000 i=1,kolik
        call kus(radka,k,t15)
        if(index(t15,'.').gt.0) call RealToFreeForm(t15)
        idl=idel(t15)
        if(i.eq.1) then
          t80=t15(1:idl)
        else
          t80=t80(1:idel(t80))//' '//t15(1:idl)
        endif
2000  continue
      radka=t80
      return
      end
      subroutine RealToFreeForm(cislo)
      implicit real*8 (a-h,o-z)
      character*(*) cislo
      idl=idel(cislo)
1000  if(cislo(idl:idl).eq.'0') then
        cislo(idl:idl)=' '
        idl=idl-1
        if(idl.gt.0) go to 1000
      endif
      if(cislo(idl:idl).eq.'.') cislo(idl:idl)=' '
      if(cislo.eq.' ') cislo='0'
      return
      end
      function scalmul(u,vg)
      implicit real*8 (a-h,o-z)
      dimension u(*),vg(*)
      scalmul=0.
      do 1000 i=1,3
        scalmul=scalmul+u(i)*vg(i)
1000  continue
      return
      end
      block data Basic
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      data dta,out/5,6/,ndim,ndimi,ndimq/3,0,9/
      data (iga(j),ipga(j),idla(j),grupaa(j),itxta(j),j=1,99) &
           /  1, 1,1,'P1      ','P1                            ', &
              2, 2,1,'P-1     ','-P1                           ', &
              3, 3,1,'P2      ','P2x                           ', &
              3, 3,1,'P2      ','P2y                           ', &
              3, 3,1,'P2      ','P2z                           ', &
              4, 3,1,'P21     ','P2xa                          ', &
              4, 3,1,'P21     ','P2yb                          ', &
              4, 3,1,'P21     ','P2zc                          ', &
              5, 3,1,'C2      ','C2x                           ', &
              5, 3,1,'A2      ','A2y                           ', &
              5, 3,1,'B2      ','B2z                           ', &
              5, 3,1,'B2      ','B2x                           ', &
              5, 3,1,'C2      ','C2y                           ', &
              5, 3,1,'A2      ','A2z                           ', &
              6, 4,1,'Pm      ','P-2x                          ', &
              6, 4,1,'Pm      ','P-2y                          ', &
              6, 4,1,'Pm      ','P-2z                          ', &
              7, 4,1,'Pc      ','P-2xc                         ', &
              7, 4,1,'Pa      ','P-2ya                         ', &
              7, 4,1,'Pb      ','P-2zb                         ', &
              7, 4,1,'Pb      ','P-2xb                         ', &
              7, 4,1,'Pc      ','P-2yc                         ', &
              7, 4,1,'Pa      ','P-2za                         ', &
              7, 4,1,'Pn      ','P-2xbc                        ', &
              7, 4,1,'Pn      ','P-2yac                        ', &
              7, 4,1,'Pn      ','P-2zab                        ', &
              8, 4,1,'Cm      ','C-2x                          ', &
              8, 4,1,'Am      ','A-2y                          ', &
              8, 4,1,'Bm      ','B-2z                          ', &
              8, 4,1,'Bm      ','B-2x                          ', &
              8, 4,1,'Cm      ','C-2y                          ', &
              8, 4,1,'Am      ','A-2z                          ', &
              9, 4,1,'Cc      ','C-2xc                         ', &
              9, 4,1,'Aa      ','A-2ya                         ', &
              9, 4,1,'Bb      ','B-2zb                         ', &
              9, 4,1,'Bb      ','B-2xb                         ', &
              9, 4,1,'Cc      ','C-2yc                         ', &
              9, 4,1,'Aa      ','A-2za                         ', &
             10, 5,2,'P2/m    ','-P2x                          ', &
             10, 5,2,'P2/m    ','-P2y                          ', &
             10, 5,2,'P2/m    ','-P2z                          ', &
             11, 5,2,'P21/m   ','-P2xa                         ', &
             11, 5,2,'P21/m   ','-P2yb                         ', &
             11, 5,2,'P21/m   ','-P2zc                         ', &
             12, 5,2,'C2/m    ','-C2x                          ', &
             12, 5,2,'A2/m    ','-A2y                          ', &
             12, 5,2,'B2/m    ','-B2z                          ', &
             12, 5,2,'B2/m    ','-B2x                          ', &
             12, 5,2,'C2/m    ','-C2y                          ', &
             12, 5,2,'A2/m    ','-A2z                          ', &
             13, 5,2,'P2/c    ','-P2xc                         ', &
             13, 5,2,'P2/a    ','-P2ya                         ', &
             13, 5,2,'P2/b    ','-P2zb                         ', &
             13, 5,2,'P2/b    ','-P2xb                         ', &
             13, 5,2,'P2/c    ','-P2yc                         ', &
             13, 5,2,'P2/a    ','-P2za                         ', &
             13, 5,2,'P2/n    ','-P2xbc                        ', &
             13, 5,2,'P2/n    ','-P2yac                        ', &
             13, 5,2,'P2/n    ','-P2zab                        ', &
             14, 5,2,'P21/c   ','-P2xca                        ', &
             14, 5,2,'P21/a   ','-P2yab                        ', &
             14, 5,2,'P21/b   ','-P2zbc                        ', &
             14, 5,2,'P21/b   ','-P2xba                        ', &
             14, 5,2,'P21/c   ','-P2ycb                        ', &
             14, 5,2,'P21/a   ','-P2zac                        ', &
             14, 5,2,'P21/n   ','-P2xabc                       ', &
             14, 5,2,'P21/n   ','-P2yabc                       ', &
             14, 5,2,'P21/n   ','-P2zabc                       ', &
             15, 5,2,'C2/c    ','-C2xc                         ', &
             15, 5,2,'A2/a    ','-A2ya                         ', &
             15, 5,2,'B2/b    ','-B2zb                         ', &
             15, 5,2,'B2/b    ','-B2xb                         ', &
             15, 5,2,'C2/c    ','-C2yc                         ', &
             15, 5,2,'A2/a    ','-A2za                         ', &
             16, 6,3,'P222    ','P2z;2x                        ', &
             17, 6,3,'P2221   ','P2zc;2x                       ', &
             17, 6,3,'P2122   ','P2xa;2y                       ', &
             17, 6,3,'P2212   ','P2yb;2z                       ', &
             18, 6,3,'P21212  ','P2z;2xab                      ', &
             18, 6,3,'P22121  ','P2x;2ybc                      ', &
             18, 6,3,'P21221  ','P2y;2zca                      ', &
             19, 6,3,'P212121 ','P2zac;2xab                    ', &
             20, 6,3,'C2221   ','C2zc;2x                       ', &
             20, 6,3,'A2122   ','A2xa;2y                       ', &
             20, 6,3,'B2212   ','B2yb;2z                       ', &
             21, 6,3,'C222    ','C2z;2x                        ', &
             21, 6,3,'A222    ','A2x;2y                        ', &
             21, 6,3,'B222    ','B2y;2z                        ', &
             22, 6,3,'F222    ','F2z;2x                        ', &
             23, 6,3,'I222    ','I2z;2x                        ', &
             24, 6,3,'I212121 ','I2zac;2xab                    ', &
             25, 7,3,'Pmm2    ','P2z;-2x                       ', &
             25, 7,3,'P2mm    ','P2x;-2y                       ', &
             25, 7,3,'Pm2m    ','P2y;-2z                       ', &
             26, 7,3,'Pmc21   ','P2zc;-2x                      ', &
             26, 7,3,'P21ma   ','P2xa;-2y                      ', &
             26, 7,3,'Pb21m   ','P2yb;-2z                      ', &
             26, 7,3,'Pcm21   ','P2zc;-2y                      ', &
             26, 7,3,'P21am   ','P2xa;-2z                      '/
!
      data (iga(j),ipga(j),idla(j),grupaa(j),itxta(j),j=100,109) &
           / 26, 7,3,'Pm21b   ','P2yb;-2x                      ', &
             27, 7,3,'Pcc2    ','P2z;-2xc                      ', &
             27, 7,3,'P2aa    ','P2x;-2ya                      ', &
             27, 7,3,'Pb2b    ','P2y;-2zb                      ', &
             28, 7,3,'Pma2    ','P2z;-2xa                      ', &
             28, 7,3,'P2mb    ','P2x;-2yb                      ', &
             28, 7,3,'Pc2m    ','P2y;-2zc                      ', &
             28, 7,3,'Pbm2    ','P2z;-2yb                      ', &
             28, 7,3,'P2cm    ','P2x;-2zc                      ', &
             28, 7,3,'Pm2a    ','P2y;-2xa                      '/
      data (iga(j),ipga(j),idla(j),grupaa(j),itxta(j),j=110,198) &
           / 29, 7,3,'Pca21   ','P2zc;-2xac                    ', &
             29, 7,3,'P21ab   ','P2xa;-2yba                    ', &
             29, 7,3,'Pc21b   ','P2yb;-2zcb                    ', &
             29, 7,3,'Pbc21   ','P2zc;-2ybc                    ', &
             29, 7,3,'P21ca   ','P2xa;-2zca                    ', &
             29, 7,3,'Pb21a   ','P2yb;-2xab                    ', &
             30, 7,3,'Pnc2    ','P2z;-2xbc                     ', &
             30, 7,3,'P2na    ','P2x;-2yca                     ', &
             30, 7,3,'Pb2n    ','P2y;-2zab                     ', &
             30, 7,3,'Pcn2    ','P2z;-2yac                     ', &
             30, 7,3,'P2an    ','P2x;-2zba                     ', &
             30, 7,3,'Pn2b    ','P2y;-2xcb                     ', &
             31, 7,3,'Pmn21   ','P2zac;-2x                     ', &
             31, 7,3,'P21mn   ','P2xba;-2y                     ', &
             31, 7,3,'Pn21m   ','P2ycb;-2z                     ', &
             31, 7,3,'Pnm21   ','P2zbc;-2y                     ', &
             31, 7,3,'P21nm   ','P2xca;-2z                     ', &
             31, 7,3,'Pm21n   ','P2yab;-2x                     ', &
             32, 7,3,'Pba2    ','P2z;-2xab                     ', &
             32, 7,3,'P2cb    ','P2x;-2ybc                     ', &
             32, 7,3,'Pc2a    ','P2y;-2zca                     ', &
             33, 7,3,'Pna21   ','P2zc;-2xn                     ', &
             33, 7,3,'P21nb   ','P2xa;-2yn                     ', &
             33, 7,3,'Pc21n   ','P2yb;-2zn                     ', &
             33, 7,3,'Pbn21   ','P2zc;-2yn                     ', &
             33, 7,3,'P21cn   ','P2xa;-2zn                     ', &
             33, 7,3,'Pn21a   ','P2yb;-2xn                     ', &
             34, 7,3,'Pnn2    ','P2z;-2xn                      ', &
             34, 7,3,'P2nn    ','P2x;-2yn                      ', &
             34, 7,3,'Pn2n    ','P2y;-2zn                      ', &
             35, 7,3,'Cmm2    ','C2z;-2x                       ', &
             35, 7,3,'A2mm    ','A2x;-2y                       ', &
             35, 7,3,'Bm2m    ','B2y;-2z                       ', &
             36, 7,3,'Cmc21   ','C2zc;-2x                      ', &
             36, 7,3,'A21ma   ','A2xa;-2y                      ', &
             36, 7,3,'Bb21m   ','B2yb;-2z                      ', &
             36, 7,3,'Ccm21   ','C2zc;-2y                      ', &
             36, 7,3,'A21am   ','A2xa;-2z                      ', &
             36, 7,3,'Bm21b   ','B2yb;-2x                      ', &
             37, 7,3,'Ccc2    ','C2z;-2xc                      ', &
             37, 7,3,'A2aa    ','A2x;-2ya                      ', &
             37, 7,3,'Bb2b    ','B2y;-2zb                      ', &
             38, 7,3,'Amm2    ','A2z;-2x                       ', &
             38, 7,3,'B2mm    ','B2x;-2y                       ', &
             38, 7,3,'Cm2m    ','C2y;-2z                       ', &
             39, 7,3,'Abm2    ','A2z;-2xb                      ', &
             39, 7,3,'B2cm    ','B2x;-2yc                      ', &
             39, 7,3,'Cm2a    ','C2y;-2za                      ', &
             39, 7,3,'Bma2    ','B2z;-2ya                      ', &
             39, 7,3,'C2mb    ','C2x;-2zb                      ', &
             39, 7,3,'Ac2m    ','A2y;-2xc                      ', &
             40, 7,3,'Ama2    ','A2z;-2xa                      ', &
             40, 7,3,'B2mb    ','B2x;-2yb                      ', &
             40, 7,3,'Cc2m    ','C2y;-2zc                      ', &
             40, 7,3,'Bbm2    ','B2z;-2yb                      ', &
             40, 7,3,'C2cm    ','C2x;-2zc                      ', &
             40, 7,3,'Am2a    ','A2y;-2xa                      ', &
             41, 7,3,'Aba2    ','A2z;-2xab                     ', &
             41, 7,3,'B2cb    ','B2x;-2ybc                     ', &
             41, 7,3,'Cc2a    ','C2y;-2zca                     ', &
             41, 7,3,'Bba2    ','B2z;-2yba                     ', &
             41, 7,3,'C2cb    ','C2x;-2zcb                     ', &
             41, 7,3,'Ac2a    ','A2y;-2xac                     ', &
             42, 7,3,'Fmm2    ','F2z;-2x                       ', &
             42, 7,3,'F2mm    ','F2x;-2y                       ', &
             42, 7,3,'Fm2m    ','F2y;-2z                       ', &
             43, 7,3,'Fdd2    ','F2z;-2xd                      ', &
             43, 7,3,'F2dd    ','F2x;-2yd                      ', &
             43, 7,3,'Fd2d    ','F2y;-2zd                      ', &
             44, 7,3,'Imm2    ','I2z;-2x                       ', &
             44, 7,3,'I2mm    ','I2x;-2y                       ', &
             44, 7,3,'Im2m    ','I2y;-2z                       ', &
             45, 7,3,'Iba2    ','I2z;-2xab                     ', &
             45, 7,3,'I2cb    ','I2x;-2ybc                     ', &
             45, 7,3,'Ic2a    ','I2y;-2zca                     ', &
             46, 7,3,'Ima2    ','I2z;-2xa                      ', &
             46, 7,3,'I2mb    ','I2x;-2yb                      ', &
             46, 7,3,'Ic2m    ','I2y;-2zc                      ', &
             46, 7,3,'Ibm2    ','I2z;-2yb                      ', &
             46, 7,3,'I2cm    ','I2x;-2zc                      ', &
             46, 7,3,'Im2a    ','I2y;-2xa                      ', &
             47, 8,3,'Pmmm    ','-P-2z;-2x                     ', &
             48, 8,3,'Pnnn    ','-P-2zab;-2xbc                 ', &
             49, 8,3,'Pccm    ','-P-2z;-2xc                    ', &
             49, 8,3,'Pmaa    ','-P-2x;-2ya                    ', &
             49, 8,3,'Pbmb    ','-P-2y;-2zb                    ', &
             50, 8,3,'Pban    ','-P-2zab;-2xb                  ', &
             50, 8,3,'Pncb    ','-P-2xbc;-2yc                  ', &
             50, 8,3,'Pcna    ','-P-2yca;-2za                  '/
      data (iga(j),ipga(j),idla(j),grupaa(j),itxta(j),j=199,297) &
           / 51, 8,3,'Pmma    ','-P-2za;-2xa                   ', &
             51, 8,3,'Pbmm    ','-P-2xb;-2yb                   ', &
             51, 8,3,'Pmcm    ','-P-2yc;-2zc                   ', &
             51, 8,3,'Pmam    ','-P-2ya;-2xa                   ', &
             51, 8,3,'Pmmb    ','-P-2zb;-2yb                   ', &
             51, 8,3,'Pcmm    ','-P-2xc;-2zc                   ', &
             52, 8,3,'Pnna    ','-P-2za;-2xbc                  ', &
             52, 8,3,'Pbnn    ','-P-2xb;-2yca                  ', &
             52, 8,3,'Pncn    ','-P-2yc;-2zab                  ', &
             52, 8,3,'Pnan    ','-P-2ya;-2xbc                  ', &
             52, 8,3,'Pnnb    ','-P-2zb;-2yca                  ', &
             52, 8,3,'Pcnn    ','-P-2xc;-2zab                  ', &
             53, 8,3,'Pmna    ','-P-2zac;-2x                   ', &
             53, 8,3,'Pbmn    ','-P-2xba;-2y                   ', &
             53, 8,3,'Pncm    ','-P-2ycb;-2z                   ', &
             53, 8,3,'Pman    ','-P-2yab;-2x                   ', &
             53, 8,3,'Pnmb    ','-P-2zbc;-2y                   ', &
             53, 8,3,'Pcnm    ','-P-2xca;-2z                   ', &
             54, 8,3,'Pcca    ','-P-2za;-2xac                  ', &
             54, 8,3,'Pbaa    ','-P-2xb;-2yba                  ', &
             54, 8,3,'Pbcb    ','-P-2yc;-2zcb                  ', &
             54, 8,3,'Pbab    ','-P-2ya;-2xab                  ', &
             54, 8,3,'Pccb    ','-P-2zb;-2ybc                  ', &
             54, 8,3,'Pcaa    ','-P-2xc;-2zca                  ', &
             55, 8,3,'Pbam    ','-P-2z;-2xab                   ', &
             55, 8,3,'Pmcb    ','-P-2x;-2ybc                   ', &
             55, 8,3,'Pcma    ','-P-2y;-2zca                   ', &
             56, 8,3,'Pccn    ','-P-2zab;-2xac                 ', &
             56, 8,3,'Pnaa    ','-P-2xbc;-2yba                 ', &
             56, 8,3,'Pbnb    ','-P-2yca;-2zcb                 ', &
             57, 8,3,'Pbcm    ','-P-2zc;-2xb                   ', &
             57, 8,3,'Pmca    ','-P-2xa;-2yc                   ', &
             57, 8,3,'Pbma    ','-P-2yb;-2za                   ', &
             57, 8,3,'Pcmb    ','-P-2yb;-2xc                   ', &
             57, 8,3,'Pcam    ','-P-2zc;-2ya                   ', &
             57, 8,3,'Pmab    ','-P-2xa;-2zb                   ', &
             58, 8,3,'Pnnm    ','-P-2z;-2xn                    ', &
             58, 8,3,'Pmnn    ','-P-2x;-2yn                    ', &
             58, 8,3,'Pnmn    ','-P-2y;-2zn                    ', &
             59, 8,3,'Pmmn    ','-P-2zab;-2xa                  ', &
             59, 8,3,'Pnmm    ','-P-2xbc;-2yb                  ', &
             59, 8,3,'Pmnm    ','-P-2yca;-2zc                  ', &
             60, 8,3,'Pbcn    ','-P-2zn;-2xab                  ', &
             60, 8,3,'Pnca    ','-P-2xn;-2ybc                  ', &
             60, 8,3,'Pbna    ','-P-2yn;-2zca                  ', &
             60, 8,3,'Pcnb    ','-P-2yn;-2xac                  ', &
             60, 8,3,'Pcan    ','-P-2zn;-2yba                  ', &
             60, 8,3,'Pnab    ','-P-2xn;-2zcb                  ', &
             61, 8,3,'Pbca    ','-P-2zac;-2xab                 ', &
             61, 8,3,'Pcab    ','-P-2yab;-2xac                 ', &
             62, 8,3,'Pnma    ','-P-2zac;-2xn                  ', &
             62, 8,3,'Pbnm    ','-P-2xba;-2yn                  ', &
             62, 8,3,'Pmcn    ','-P-2ycb;-2zn                  ', &
             62, 8,3,'Pnam    ','-P-2yab;-2xn                  ', &
             62, 8,3,'Pmnb    ','-P-2zbc;-2yn                  ', &
             62, 8,3,'Pcmn    ','-P-2xca;-2zn                  ', &
             63, 8,3,'Cmcm    ','-C-2zc;-2x                    ', &
             63, 8,3,'Amma    ','-A-2xa;-2y                    ', &
             63, 8,3,'Bbmm    ','-B-2yb;-2z                    ', &
             63, 8,3,'Bmmb    ','-B-2yb;-2x                    ', &
             63, 8,3,'Ccmm    ','-C-2zc;-2y                    ', &
             63, 8,3,'Amam    ','-A-2xa;-2z                    ', &
             64, 8,3,'Cmca    ','-C-2zac;-2x                   ', &
             64, 8,3,'Abma    ','-A-2xba;-2y                   ', &
             64, 8,3,'Bbcm    ','-B-2ycb;-2z                   ', &
             64, 8,3,'Bmab    ','-B-2yab;-2x                   ', &
             64, 8,3,'Ccmb    ','-C-2zbc;-2y                   ', &
             64, 8,3,'Acam    ','-A-2xca;-2z                   ', &
             65, 8,3,'Cmmm    ','-C-2z;-2x                     ', &
             65, 8,3,'Ammm    ','-A-2x;-2y                     ', &
             65, 8,3,'Bmmm    ','-B-2y;-2z                     ', &
             66, 8,3,'Cccm    ','-C-2z;-2xc                    ', &
             66, 8,3,'Amaa    ','-A-2x;-2ya                    ', &
             66, 8,3,'Bbmb    ','-B-2y;-2zb                    ', &
             67, 8,3,'Cmma    ','-C-2za;-2x                    ', &
             67, 8,3,'Abmm    ','-A-2xb;-2y                    ', &
             67, 8,3,'Bmcm    ','-B-2yc;-2z                    ', &
             67, 8,3,'Bmam    ','-B-2ya;-2x                    ', &
             67, 8,3,'Cmmb    ','-C-2zb;-2y                    ', &
             67, 8,3,'Acmm    ','-A-2xc;-2z                    ', &
             68, 8,3,'Ccca    ','-C-2za;-2xac                  ', &
             68, 8,3,'Abaa    ','-A-2xb;-2yba                  ', &
             68, 8,3,'Bbcb    ','-B-2yc;-2zcb                  ', &
             68, 8,3,'Bbab    ','-B-2ya;-2xab                  ', &
             68, 8,3,'Cccb    ','-C-2zb;-2ybc                  ', &
             68, 8,3,'Acaa    ','-A-2xc;-2zca                  ', &
             69, 8,3,'Fmmm    ','-F-2z;-2x                     ', &
             70, 8,3,'Fddd    ','-F-2zuv;-2xvw                 ', &
             71, 8,3,'Immm    ','-I-2z;-2x                     ', &
             72, 8,3,'Ibam    ','-I-2z;-2xab                   ', &
             72, 8,3,'Imcb    ','-I-2x;-2ybc                   ', &
             72, 8,3,'Icma    ','-I-2y;-2zca                   ', &
             73, 8,3,'Ibca    ','-I-2zac;-2xab                 ', &
             73, 8,3,'Icab    ','-I-2yab;-2xac                 ', &
             74, 8,3,'Imma    ','-I-2zac;-2x                   ', &
             74, 8,3,'Ibmm    ','-I-2xba;-2y                   ', &
             74, 8,3,'Imcm    ','-I-2ycb;-2z                   ', &
             74, 8,3,'Imam    ','-I-2yab;-2x                   ', &
             74, 8,3,'Immb    ','-I-2zbc;-2y                   '/
      data (iga(j),ipga(j),idla(j),grupaa(j),itxta(j),j=298,307) &
           / 74, 8,3,'Icmm    ','-I-2xca;-2z                   ', &
             75, 9,1,'P4      ','P4                            ', &
             76, 9,1,'P41     ','P41                           ', &
             77, 9,1,'P42     ','P4c                           ', &
             78, 9,1,'P43     ','P43                           ', &
             79, 9,1,'I4      ','I4                            ', &
             80, 9,1,'I41     ','I41b                          ', &
             81,10,1,'P-4     ','P-4                           ', &
             82,10,1,'I-4     ','I-4                           ', &
             83,11,2,'P4/m    ','-P4                           '/
      data (iga(j),ipga(j),idla(j),grupaa(j),itxta(j),j=308,396) &
           / 84,11,2,'P42/m   ','-P4c                          ', &
             85,11,2,'P4/n    ','-P4a                          ', &
             86,11,2,'P42/n   ','-P4bc                         ', &
             87,11,2,'I4/m    ','-I4                           ', &
             88,11,2,'I41/a   ','-I4ad                         ', &
             89,12,3,'P422    ','P4;2                          ', &
             90,12,3,'P4212   ','P4ab;2ab                      ', &
             91,12,3,'P4122   ','P41;2c                        ', &
             92,12,3,'P41212  ','P43n;2nw                      ', &
             93,12,3,'P4222   ','P4c;2                         ', &
             94,12,3,'P42212  ','P4n;2n                        ', &
             95,12,3,'P4322   ','P43;2c                        ', &
             96,12,3,'P43212  ','P41n;2abw                     ', &
             97,12,3,'I422    ','I4;2                          ', &
             98,12,3,'I4122   ','I41b;2bw                      ', &
             99,13,3,'P4mm    ','P4;-2                         ', &
            100,13,3,'P4bm    ','P4;-2ab                       ', &
            101,13,3,'P42cm   ','P4c;-2c                       ', &
            102,13,3,'P42nm   ','P4n;-2n                       ', &
            103,13,3,'P4cc    ','P4;-2c                        ', &
            104,13,3,'P4nc    ','P4;-2n                        ', &
            105,13,3,'P42mc   ','P4c;-2                        ', &
            106,13,3,'P42bc   ','P4c;-2ab                      ', &
            107,13,3,'I4mm    ','I4;-2                         ', &
            108,13,3,'I4cm    ','I4;-2ab                       ', &
            109,13,3,'I41md   ','I41b;-2                       ', &
            110,13,3,'I41cd   ','I41b;-2c                      ', &
            111,14,3,'P-42m   ','P-4;2                         ', &
            112,14,3,'P-42c   ','P-4;2c                        ', &
            113,14,3,'P-421m  ','P-4;2ab                       ', &
            114,14,3,'P-421c  ','P-4;2n                        ', &
            115,14,3,'P-4m2   ','P-4;-2                        ', &
            116,14,3,'P-4c2   ','P-4;-2c                       ', &
            117,14,3,'P-4b2   ','P-4;-2ab                      ', &
            118,14,3,'P-4n2   ','P-4;-2n                       ', &
            119,14,3,'I-4m2   ','I-4;-2                        ', &
            120,14,3,'I-4c2   ','I-4;-2c                       ', &
            121,14,3,'I-42m   ','I-4;2                         ', &
            122,14,3,'I-42d   ','I-4;2bw                       ', &
            123,15,4,'P4/mmm  ','-P4;-2                        ', &
            124,15,4,'P4/mcc  ','-P4;-2c                       ', &
            125,15,4,'P4/nbm  ','-P4a;-2b                      ', &
            126,15,4,'P4/nnc  ','-P4a;-2bc                     ', &
            127,15,4,'P4/mbm  ','-P4;-2ab                      ', &
            128,15,4,'P4/mnc  ','-P4;-2n                       ', &
            129,15,4,'P4/nmm  ','-P4a;-2a                      ', &
            130,15,4,'P4/ncc  ','-P4a;-2ac                     ', &
            131,15,4,'P42/mmc ','-P4c;-2                       ', &
            132,15,4,'P42/mcm ','-P4c;-2c                      ', &
            133,15,4,'P42/nbc ','-P4ac;-2b                     ', &
            134,15,4,'P42/nnm ','-P4ac;-2bc                    ', &
            135,15,4,'P42/mbc ','-P4c;-2ab                     ', &
            136,15,4,'P42/mnm ','-P4n;-2n                      ', &
            137,15,4,'P42/nmc ','-P4ac;-2a                     ', &
            138,15,4,'P42/ncm ','-P4ac;-2ac                    ', &
            139,15,4,'I4/mmm  ','-I4;-2                        ', &
            140,15,4,'I4/mcm  ','-I4;-2c                       ', &
            141,15,4,'I41/amd ','-I4bd;-2                      ', &
            142,15,4,'I41/acd ','-I4bd;-2c                     ', &
            143,16,1,'P3      ','P3                            ', &
            144,16,1,'P31     ','P31                           ', &
            145,16,1,'P32     ','P32                           ', &
            146,16,1,'R3      ','R3                            ', &
            147,17,1,'P-3     ','-P3                           ', &
            148,17,1,'R-3     ','-R3                           ', &
            149,18,3,'P312    ','P3;2                          ', &
            150,18,3,'P321    ','P3;2"                         ', &
            151,18,3,'P3112   ','P31;2#0,0,1/3                 ', &
            152,18,3,'P3121   ','P31;2"                        ', &
            153,18,3,'P3212   ','P32;2#0,0,1/6                 ', &
            154,18,3,'P3221   ','P32;2"                        ', &
            155,18,2,'R32     ','R3;2"                         ', &
            156,19,3,'P3m1    ','P3;-2"                        ', &
            157,19,3,'P31m    ','P3;-2                         ', &
            158,19,3,'P3c1    ','P3;-2"c                       ', &
            159,19,3,'P31c    ','P3;-2c                        ', &
            160,19,2,'R3m     ','R3;-2"                        ', &
            161,19,2,'R3c     ','R3;-2"c                       ', &
            162,20,3,'P-31m   ','-P3;-2                        ', &
            163,20,3,'P-31c   ','-P3;-2c                       ', &
            164,20,3,'P-3m1   ','-P3;-2"                       ', &
            165,20,3,'P-3c1   ','-P3;-2"c                      ', &
            166,20,2,'R-3m    ','-R3;-2"                       ', &
            167,20,2,'R-3c    ','-R3;-2"c                      ', &
            168,21,1,'P6      ','P6                            ', &
            169,21,1,'P61     ','P61                           ', &
            170,21,1,'P65     ','P65                           ', &
            171,21,1,'P62     ','P62                           ', &
            172,21,1,'P64     ','P64                           '/
      data (iga(j),ipga(j),idla(j),grupaa(j),itxta(j),j=397,454) &
           /173,21,1,'P63     ','P6c                           ', &
            174,22,1,'P-6     ','P-6                           ', &
            175,23,2,'P6/m    ','-P6                           ', &
            176,23,2,'P63/m   ','-P6c                          ', &
            177,24,3,'P622    ','P6;2                          ', &
            178,24,3,'P6122   ','P61;2#0,0,-1/12               ', &
            179,24,3,'P6522   ','P65;2#0,0,1/12                ', &
            180,24,3,'P6222   ','P62;2#0,0,1/3                 ', &
            181,24,3,'P6422   ','P64;2#0,0,1/6                 ', &
            182,24,3,'P6322   ','P6c;2#0,0,1/4                 ', &
            183,25,3,'P6mm    ','P6;-2                         ', &
            184,25,3,'P6cc    ','P6;-2c                        ', &
            185,25,3,'P63cm   ','P6c;-2                        ', &
            186,25,3,'P63mc   ','P6c;-2c                       ', &
            187,26,3,'P-6m2   ','P-6;2                         ', &
            188,26,3,'P-6c2   ','P-6c;2                        ', &
            189,26,3,'P-62m   ','P-6;-2                        ', &
            190,26,3,'P-62c   ','P-6c;-2c                      ', &
            191,27,4,'P6/mmm  ','-P6;-2                        ', &
            192,27,4,'P6/mcc  ','-P6;-2c                       ', &
            193,27,4,'P63/mcm ','-P6c;-2                       ', &
            194,27,4,'P63/mmc ','-P6c;-2c                      ', &
            195,28,2,'P23     ','P2;2;3                        ', &
            196,28,2,'F23     ','F2;2;3                        ', &
            197,28,2,'I23     ','I2;2;3                        ', &
            198,28,2,'P213    ','P2ac;2ab;3                    ', &
            199,28,2,'I213    ','I2ac;2ab;3                    ', &
            200,29,2,'Pm-3    ','-P2;2;3                       ', &
            201,29,2,'Pn-3    ','-P2ab;2bc;3                   ', &
            202,29,2,'Fm-3    ','-F2;2;3                       ', &
            203,29,2,'Fd-3    ','-F2uv;2vw;3                   ', &
            204,29,2,'Im-3    ','-I2;2;3                       ', &
            205,29,2,'Pa-3    ','-P2ac;2ab;3                   ', &
            206,29,2,'Ia-3    ','-I2ac;2ab;3                   ', &
            207,30,3,'P432    ','P4;2;3                        ', &
            208,30,3,'P4232   ','P4n;2;3                       ', &
            209,30,3,'F432    ','F4;2;3                        ', &
            210,30,3,'F4132   ','F4d;2;3                       ', &
            211,30,3,'I432    ','I4;2;3                        ', &
            212,30,3,'P4332   ','P4bdn;2ab;3                   ', &
            213,30,3,'P4132   ','P4bd;2ab;3                    ', &
            214,30,3,'I4132   ','I4bd;2ab;3                    ', &
            215,31,3,'P-43m   ','P-4;2;3                       ', &
            216,31,3,'F-43m   ','F-4;2;3                       ', &
            217,31,3,'I-43m   ','I-4;2;3                       ', &
            218,31,3,'P-43n   ','P-4n;2;3                      ', &
            219,31,3,'F-43c   ','F-4c;2;3                      ', &
            220,31,3,'I-43d   ','I-4bd;2ab;3                   ', &
            221,32,3,'Pm-3m   ','-P4;2;3                       ', &
            222,32,3,'Pn-3n   ','-P4a;2bc;3                    ', &
            223,32,3,'Pm-3n   ','-P4n;2;3                      ', &
            224,32,3,'Pn-3m   ','-P4bc;2bc;3                   ', &
            225,32,3,'Fm-3m   ','-F4;2;3                       ', &
            226,32,3,'Fm-3c   ','-F4n;2;3                      ', &
            227,32,3,'Fd-3m   ','-F4vw;2vw;3                   ', &
            228,32,3,'Fd-3c   ','-F4ud;2vw;3                   ', &
            229,32,3,'Im-3m   ','-I4;2;3                       ', &
            230,32,3,'Ia-3d   ','-I4bd;2ab;3                   '/
      data smbx/'x','y','z'/
      data smbc/'PABCIRF'/,FormA80,FormA1/'(a80)','(333a1)'/
      end
