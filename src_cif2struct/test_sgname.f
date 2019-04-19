 subroutine test_sgname(sgname1,sopexist)

   use struct, only: sgname,sgnum,lattyp
   character*32  sgname1
   logical       sopexist

   call ChangeToTitle(sgname1)
   call Zhusti(sgname1)
   sgname=sgname1
 
   call getlattype(sopexist)

   sgname1=sgname
   call getsgnum


 end subroutine test_sgname

