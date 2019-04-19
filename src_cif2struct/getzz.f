 subroutine getzz

   use struct, only : aname,nat,r0,zz
   character elem*2,text*80
   real*8 elmzz
   integer i

   do i=1,nat

1    continue
     elem='  '
     elem=trim(aname(i)(1:2))
     elmzz=-100.0d0
      
     if (elem.eq.'Ac') elmzz= 89. 
     if (elem.eq.'Ag') elmzz= 47. 
     if (elem.eq.'Al') elmzz= 13. 
     if (elem.eq.'Am') elmzz= 95. 
     if (elem.eq.'Ar') elmzz= 18. 
     if (elem.eq.'As') elmzz= 33. 
     if (elem.eq.'At') elmzz= 85. 
     if (elem.eq.'Au') elmzz= 79. 
     if (elem.eq.'B') elmzz=  5. 
     if (elem.eq.'Ba') elmzz= 56. 
     if (elem.eq.'Be') elmzz=  4. 
     if (elem.eq.'Bi') elmzz= 83. 
     if (elem.eq.'Bk') elmzz= 97. 
     if (elem.eq.'Br') elmzz= 35. 
     if (elem.eq.'C') elmzz=  6. 
     if (elem.eq.'Ca') elmzz= 20. 
     if (elem.eq.'Cd') elmzz= 48. 
     if (elem.eq.'Ce') elmzz= 58. 
     if (elem.eq.'Cf') elmzz= 98. 
     if (elem.eq.'Cl') elmzz= 17. 
     if (elem.eq.'Cm') elmzz= 96. 
     if (elem.eq.'Co') elmzz= 27. 
     if (elem.eq.'Cr') elmzz= 24. 
     if (elem.eq.'Cs') elmzz= 55. 
     if (elem.eq.'Cu') elmzz= 29. 
     if (elem.eq.'Dy') elmzz= 66. 
     if (elem.eq.'Er') elmzz= 68. 
     if (elem.eq.'Es') elmzz= 99. 
     if (elem.eq.'Eu') elmzz= 63. 
     if (elem.eq.'F') elmzz=  9. 
     if (elem.eq.'Fe') elmzz= 26. 
     if (elem.eq.'Fm') elmzz=100. 
     if (elem.eq.'Fr') elmzz= 87. 
     if (elem.eq.'Ga') elmzz= 31. 
     if (elem.eq.'Gd') elmzz= 64. 
     if (elem.eq.'Ge') elmzz= 32. 
     if (elem.eq.'H'.OR.elem.eq.'D') elmzz=  1. 
     if (elem.eq.'He') elmzz=  2. 
     if (elem.eq.'Hf') elmzz= 72. 
     if (elem.eq.'Hg') elmzz= 80. 
     if (elem.eq.'Ho') elmzz= 67. 
     if (elem.eq.'I') elmzz= 53. 
     if (elem.eq.'In') elmzz= 49. 
     if (elem.eq.'Ir') elmzz= 77. 
     if (elem.eq.'K') elmzz= 19. 
     if (elem.eq.'Kr') elmzz= 36. 
     if (elem.eq.'La') elmzz= 57. 
     if (elem.eq.'Li') elmzz=  3. 
     if (elem.eq.'Lr') elmzz=103. 
     if (elem.eq.'Lu') elmzz= 71. 
     if (elem.eq.'Md') elmzz=101. 
     if (elem.eq.'Mg') elmzz= 12. 
     if (elem.eq.'Mn') elmzz= 25. 
     if (elem.eq.'Mo') elmzz= 42. 
     if (elem.eq.'N') elmzz=  7. 
     if (elem.eq.'Na') elmzz= 11. 
     if (elem.eq.'Nb') elmzz= 41. 
     if (elem.eq.'Nd') elmzz= 60. 
     if (elem.eq.'Ne') elmzz= 10. 
     if (elem.eq.'Ni') elmzz= 28. 
     if (elem.eq.'No') elmzz=102. 
     if (elem.eq.'Np') elmzz= 93. 
     if (elem.eq.'O') elmzz=  8. 
     if (elem.eq.'Os') elmzz= 76. 
     if (elem.eq.'P') elmzz= 15. 
     if (elem.eq.'Pa') elmzz= 91. 
     if (elem.eq.'Pb') elmzz= 82. 
     if (elem.eq.'Pd') elmzz= 46. 
     if (elem.eq.'Pm') elmzz= 61. 
     if (elem.eq.'Po') elmzz= 84. 
     if (elem.eq.'Pr') elmzz= 59. 
     if (elem.eq.'Pt') elmzz= 78. 
     if (elem.eq.'Pu') elmzz= 94. 
     if (elem.eq.'Ra') elmzz= 88. 
     if (elem.eq.'Rb') elmzz= 37. 
     if (elem.eq.'Re') elmzz= 75. 
     if (elem.eq.'Rh') elmzz= 45. 
     if (elem.eq.'Rn') elmzz= 86. 
     if (elem.eq.'Ru') elmzz= 44. 
     if (elem.eq.'S') elmzz= 16. 
     if (elem.eq.'Sb') elmzz= 51. 
     if (elem.eq.'Sc') elmzz= 21. 
     if (elem.eq.'Se') elmzz= 34. 
     if (elem.eq.'Si') elmzz= 14. 
     if (elem.eq.'Sm') elmzz= 62. 
     if (elem.eq.'Sn') elmzz= 50. 
     if (elem.eq.'Sr') elmzz= 38. 
     if (elem.eq.'Ta') elmzz= 73. 
     if (elem.eq.'Tb') elmzz= 65. 
     if (elem.eq.'Tc') elmzz= 43. 
     if (elem.eq.'Te') elmzz= 52. 
     if (elem.eq.'Th') elmzz= 90. 
     if (elem.eq.'Ti') elmzz= 22. 
     if (elem.eq.'Tl') elmzz= 81. 
     if (elem.eq.'Tm') elmzz= 69. 
     if (elem.eq.'U') elmzz= 92. 
     if (elem.eq.'V') elmzz= 23. 
     if (elem.eq.'W') elmzz= 74. 
     if (elem.eq.'Xe') elmzz= 54. 
     if (elem.eq.'Y') elmzz= 39. 
     if (elem.eq.'Yb') elmzz= 70. 
     if (elem.eq.'Zn') elmzz= 30. 
     if (elem.eq.'Zr') elmzz= 40. 

     if (elmzz.lt.0.0d0) then

        text='unknown element: '//trim(elem)//', type in correct name'
        call ReadText(trim(text),aname(i),' ')
        goto 1
     endif
     zz(i)=elmzz

     r0(i)=5.0d-6
     if (zz(i).le.71) r0(i)=1.0d-5
     if (zz(i).le.36) r0(i)=5.0d-5
     if (zz(i).le.18) r0(i)=1.0d-4
     
  enddo

end subroutine getzz
