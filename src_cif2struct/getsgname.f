
 subroutine getsgname(sgnum,sgname)

   integer sgnum
   character sgname*(*)

   sgname='        '
   if (sgnum.eq.  1) sgname='P1      '
   if (sgnum.eq.  2) sgname='P-1     '
   if (sgnum.eq.  3) sgname='P2      '
   if (sgnum.eq.  3) sgname='P2      '
   if (sgnum.eq.  3) sgname='P2      '
   if (sgnum.eq.  4) sgname='P21     '
   if (sgnum.eq.  4) sgname='P21     '
   if (sgnum.eq.  4) sgname='P21     '
   if (sgnum.eq.  5) sgname='C2      '
   if (sgnum.eq.  5) sgname='A2      '
   if (sgnum.eq.  5) sgname='B2      '
   if (sgnum.eq.  5) sgname='B2      '
   if (sgnum.eq.  5) sgname='C2      '
   if (sgnum.eq.  5) sgname='A2      '
   if (sgnum.eq.  6) sgname='Pm      '
   if (sgnum.eq.  6) sgname='Pm      '
   if (sgnum.eq.  6) sgname='Pm      '
   if (sgnum.eq.  7) sgname='Pc      '
   if (sgnum.eq.  7) sgname='Pa      '
   if (sgnum.eq.  7) sgname='Pb      '
   if (sgnum.eq.  7) sgname='Pb      '
   if (sgnum.eq.  7) sgname='Pc      '
   if (sgnum.eq.  7) sgname='Pa      '
   if (sgnum.eq.  7) sgname='Pn      '
   if (sgnum.eq.  7) sgname='Pn      '
   if (sgnum.eq.  7) sgname='Pn      '
   if (sgnum.eq.  8) sgname='Cm      '
   if (sgnum.eq.  8) sgname='Am      '
   if (sgnum.eq.  8) sgname='Bm      '
   if (sgnum.eq.  8) sgname='Bm      '
   if (sgnum.eq.  8) sgname='Cm      '
   if (sgnum.eq.  8) sgname='Am      '
   if (sgnum.eq.  9) sgname='Cc      '
   if (sgnum.eq.  9) sgname='Aa      '
   if (sgnum.eq.  9) sgname='Bb      '
   if (sgnum.eq.  9) sgname='Bb      '
   if (sgnum.eq.  9) sgname='Cc      '
   if (sgnum.eq.  9) sgname='Aa      '
   if (sgnum.eq. 10) sgname='P2/m    '
   if (sgnum.eq. 10) sgname='P2/m    '
   if (sgnum.eq. 10) sgname='P2/m    '
   if (sgnum.eq. 11) sgname='P21/m   '
   if (sgnum.eq. 11) sgname='P21/m   '
   if (sgnum.eq. 11) sgname='P21/m   '
   if (sgnum.eq. 12) sgname='C2/m    '
   if (sgnum.eq. 12) sgname='A2/m    '
   if (sgnum.eq. 12) sgname='B2/m    '
   if (sgnum.eq. 12) sgname='B2/m    '
   if (sgnum.eq. 12) sgname='C2/m    '
   if (sgnum.eq. 12) sgname='A2/m    '
   if (sgnum.eq. 13) sgname='P2/c    '
   if (sgnum.eq. 13) sgname='P2/a    '
   if (sgnum.eq. 13) sgname='P2/b    '
   if (sgnum.eq. 13) sgname='P2/b    '
   if (sgnum.eq. 13) sgname='P2/c    '
   if (sgnum.eq. 13) sgname='P2/a    '
   if (sgnum.eq. 13) sgname='P2/n    '
   if (sgnum.eq. 13) sgname='P2/n    '
   if (sgnum.eq. 13) sgname='P2/n    '
   if (sgnum.eq. 14) sgname='P21/c   '
   if (sgnum.eq. 14) sgname='P21/a   '
   if (sgnum.eq. 14) sgname='P21/b   '
   if (sgnum.eq. 14) sgname='P21/b   '
   if (sgnum.eq. 14) sgname='P21/c   '
   if (sgnum.eq. 14) sgname='P21/a   '
   if (sgnum.eq. 14) sgname='P21/n   '
   if (sgnum.eq. 14) sgname='P21/n   '
   if (sgnum.eq. 14) sgname='P21/n   '
   if (sgnum.eq. 15) sgname='C2/c    '
   if (sgnum.eq. 15) sgname='A2/a    '
   if (sgnum.eq. 15) sgname='B2/b    '
   if (sgnum.eq. 15) sgname='B2/b    '
   if (sgnum.eq. 15) sgname='C2/c    '
   if (sgnum.eq. 15) sgname='A2/a    '
   if (sgnum.eq. 16) sgname='P222    '
   if (sgnum.eq. 17) sgname='P2221   '
   if (sgnum.eq. 17) sgname='P2122   '
   if (sgnum.eq. 17) sgname='P2212   '
   if (sgnum.eq. 18) sgname='P21212  '
   if (sgnum.eq. 18) sgname='P22121  '
   if (sgnum.eq. 18) sgname='P21221  '
   if (sgnum.eq. 19) sgname='P212121 '
   if (sgnum.eq. 20) sgname='C2221   '
   if (sgnum.eq. 20) sgname='A2122   '
   if (sgnum.eq. 20) sgname='B2212   '
   if (sgnum.eq. 21) sgname='C222    '
   if (sgnum.eq. 21) sgname='A222    '
   if (sgnum.eq. 21) sgname='B222    '
   if (sgnum.eq. 22) sgname='F222    '
   if (sgnum.eq. 23) sgname='I222    '
   if (sgnum.eq. 24) sgname='I212121 '
   if (sgnum.eq. 25) sgname='Pmm2    '
   if (sgnum.eq. 25) sgname='P2mm    '
   if (sgnum.eq. 25) sgname='Pm2m    '
   if (sgnum.eq. 26) sgname='Pmc21   '
   if (sgnum.eq. 26) sgname='P21ma   '
   if (sgnum.eq. 26) sgname='Pb21m   '
   if (sgnum.eq. 26) sgname='Pcm21   '
   if (sgnum.eq. 26) sgname='P21am   '
   if (sgnum.eq. 26) sgname='Pm21b   '
   if (sgnum.eq. 27) sgname='Pcc2    '
   if (sgnum.eq. 27) sgname='P2aa    '
   if (sgnum.eq. 27) sgname='Pb2b    '
   if (sgnum.eq. 28) sgname='Pma2    '
   if (sgnum.eq. 28) sgname='P2mb    '
   if (sgnum.eq. 28) sgname='Pc2m    '
   if (sgnum.eq. 28) sgname='Pbm2    '
   if (sgnum.eq. 28) sgname='P2cm    '
   if (sgnum.eq. 28) sgname='Pm2a    '
   if (sgnum.eq. 29) sgname='Pca21   '
   if (sgnum.eq. 29) sgname='P21ab   '
   if (sgnum.eq. 29) sgname='Pc21b   '
   if (sgnum.eq. 29) sgname='Pbc21   '
   if (sgnum.eq. 29) sgname='P21ca   '
   if (sgnum.eq. 29) sgname='Pb21a   '
   if (sgnum.eq. 30) sgname='Pnc2    '
   if (sgnum.eq. 30) sgname='P2na    '
   if (sgnum.eq. 30) sgname='Pb2n    '
   if (sgnum.eq. 30) sgname='Pcn2    '
   if (sgnum.eq. 30) sgname='P2an    '
   if (sgnum.eq. 30) sgname='Pn2b    '
   if (sgnum.eq. 31) sgname='Pmn21   '
   if (sgnum.eq. 31) sgname='P21mn   '
   if (sgnum.eq. 31) sgname='Pn21m   '
   if (sgnum.eq. 31) sgname='Pnm21   '
   if (sgnum.eq. 31) sgname='P21nm   '
   if (sgnum.eq. 31) sgname='Pm21n   '
   if (sgnum.eq. 32) sgname='Pba2    '
   if (sgnum.eq. 32) sgname='P2cb    '
   if (sgnum.eq. 32) sgname='Pc2a    '
   if (sgnum.eq. 33) sgname='Pna21   '
   if (sgnum.eq. 33) sgname='P21nb   '
   if (sgnum.eq. 33) sgname='Pc21n   '
   if (sgnum.eq. 33) sgname='Pbn21   '
   if (sgnum.eq. 33) sgname='P21cn   '
   if (sgnum.eq. 33) sgname='Pn21a   '
   if (sgnum.eq. 34) sgname='Pnn2    '
   if (sgnum.eq. 34) sgname='P2nn    '
   if (sgnum.eq. 34) sgname='Pn2n    '
   if (sgnum.eq. 35) sgname='Cmm2    '
   if (sgnum.eq. 35) sgname='A2mm    '
   if (sgnum.eq. 35) sgname='Bm2m    '
   if (sgnum.eq. 36) sgname='Cmc21   '
   if (sgnum.eq. 36) sgname='A21ma   '
   if (sgnum.eq. 36) sgname='Bb21m   '
   if (sgnum.eq. 36) sgname='Ccm21   '
   if (sgnum.eq. 36) sgname='A21am   '
   if (sgnum.eq. 36) sgname='Bm21b   '
   if (sgnum.eq. 37) sgname='Ccc2    '
   if (sgnum.eq. 37) sgname='A2aa    '
   if (sgnum.eq. 37) sgname='Bb2b    '
   if (sgnum.eq. 38) sgname='Amm2    '
   if (sgnum.eq. 38) sgname='B2mm    '
   if (sgnum.eq. 38) sgname='Cm2m    '
   if (sgnum.eq. 39) sgname='Abm2    '
   if (sgnum.eq. 39) sgname='B2cm    '
   if (sgnum.eq. 39) sgname='Cm2a    '
   if (sgnum.eq. 39) sgname='Bma2    '
   if (sgnum.eq. 39) sgname='C2mb    '
   if (sgnum.eq. 39) sgname='Ac2m    '
   if (sgnum.eq. 40) sgname='Ama2    '
   if (sgnum.eq. 40) sgname='B2mb    '
   if (sgnum.eq. 40) sgname='Cc2m    '
   if (sgnum.eq. 40) sgname='Bbm2    '
   if (sgnum.eq. 40) sgname='C2cm    '
   if (sgnum.eq. 40) sgname='Am2a    '
   if (sgnum.eq. 41) sgname='Aba2    '
   if (sgnum.eq. 41) sgname='B2cb    '
   if (sgnum.eq. 41) sgname='Cc2a    '
   if (sgnum.eq. 41) sgname='Bba2    '
   if (sgnum.eq. 41) sgname='C2cb    '
   if (sgnum.eq. 41) sgname='Ac2a    '
   if (sgnum.eq. 42) sgname='Fmm2    '
   if (sgnum.eq. 42) sgname='F2mm    '
   if (sgnum.eq. 42) sgname='Fm2m    '
   if (sgnum.eq. 43) sgname='Fdd2    '
   if (sgnum.eq. 43) sgname='F2dd    '
   if (sgnum.eq. 43) sgname='Fd2d    '
   if (sgnum.eq. 44) sgname='Imm2    '
   if (sgnum.eq. 44) sgname='I2mm    '
   if (sgnum.eq. 44) sgname='Im2m    '
   if (sgnum.eq. 45) sgname='Iba2    '
   if (sgnum.eq. 45) sgname='I2cb    '
   if (sgnum.eq. 45) sgname='Ic2a    '
   if (sgnum.eq. 46) sgname='Ima2    '
   if (sgnum.eq. 46) sgname='I2mb    '
   if (sgnum.eq. 46) sgname='Ic2m    '
   if (sgnum.eq. 46) sgname='Ibm2    '
   if (sgnum.eq. 46) sgname='I2cm    '
   if (sgnum.eq. 46) sgname='Im2a    '
   if (sgnum.eq. 47) sgname='Pmmm    '
   if (sgnum.eq. 48) sgname='Pnnn    '
   if (sgnum.eq. 49) sgname='Pccm    '
   if (sgnum.eq. 49) sgname='Pmaa    '
   if (sgnum.eq. 49) sgname='Pbmb    '
   if (sgnum.eq. 50) sgname='Pban    '
   if (sgnum.eq. 50) sgname='Pncb    '
   if (sgnum.eq. 50) sgname='Pcna    '
   if (sgnum.eq. 51) sgname='Pmma    '
   if (sgnum.eq. 51) sgname='Pbmm    '
   if (sgnum.eq. 51) sgname='Pmcm    '
   if (sgnum.eq. 51) sgname='Pmam    '
   if (sgnum.eq. 51) sgname='Pmmb    '
   if (sgnum.eq. 51) sgname='Pcmm    '
   if (sgnum.eq. 52) sgname='Pnna    '
   if (sgnum.eq. 52) sgname='Pbnn    '
   if (sgnum.eq. 52) sgname='Pncn    '
   if (sgnum.eq. 52) sgname='Pnan    '
   if (sgnum.eq. 52) sgname='Pnnb    '
   if (sgnum.eq. 52) sgname='Pcnn    '
   if (sgnum.eq. 53) sgname='Pmna    '
   if (sgnum.eq. 53) sgname='Pbmn    '
   if (sgnum.eq. 53) sgname='Pncm    '
   if (sgnum.eq. 53) sgname='Pman    '
   if (sgnum.eq. 53) sgname='Pnmb    '
   if (sgnum.eq. 53) sgname='Pcnm    '
   if (sgnum.eq. 54) sgname='Pcca    '
   if (sgnum.eq. 54) sgname='Pbaa    '
   if (sgnum.eq. 54) sgname='Pbcb    '
   if (sgnum.eq. 54) sgname='Pbab    '
   if (sgnum.eq. 54) sgname='Pccb    '
   if (sgnum.eq. 54) sgname='Pcaa    '
   if (sgnum.eq. 55) sgname='Pbam    '
   if (sgnum.eq. 55) sgname='Pmcb    '
   if (sgnum.eq. 55) sgname='Pcma    '
   if (sgnum.eq. 56) sgname='Pccn    '
   if (sgnum.eq. 56) sgname='Pnaa    '
   if (sgnum.eq. 56) sgname='Pbnb    '
   if (sgnum.eq. 57) sgname='Pbcm    '
   if (sgnum.eq. 57) sgname='Pmca    '
   if (sgnum.eq. 57) sgname='Pbma    '
   if (sgnum.eq. 57) sgname='Pcmb    '
   if (sgnum.eq. 57) sgname='Pcam    '
   if (sgnum.eq. 57) sgname='Pmab    '
   if (sgnum.eq. 58) sgname='Pnnm    '
   if (sgnum.eq. 58) sgname='Pmnn    '
   if (sgnum.eq. 58) sgname='Pnmn    '
   if (sgnum.eq. 59) sgname='Pmmn    '
   if (sgnum.eq. 59) sgname='Pnmm    '
   if (sgnum.eq. 59) sgname='Pmnm    '
   if (sgnum.eq. 60) sgname='Pbcn    '
   if (sgnum.eq. 60) sgname='Pnca    '
   if (sgnum.eq. 60) sgname='Pbna    '
   if (sgnum.eq. 60) sgname='Pcnb    '
   if (sgnum.eq. 60) sgname='Pcan    '
   if (sgnum.eq. 60) sgname='Pnab    '
   if (sgnum.eq. 61) sgname='Pbca    '
   if (sgnum.eq. 61) sgname='Pcab    '
   if (sgnum.eq. 62) sgname='Pnma    '
   if (sgnum.eq. 62) sgname='Pbnm    '
   if (sgnum.eq. 62) sgname='Pmcn    '
   if (sgnum.eq. 62) sgname='Pnam    '
   if (sgnum.eq. 62) sgname='Pmnb    '
   if (sgnum.eq. 62) sgname='Pcmn    '
   if (sgnum.eq. 63) sgname='Cmcm    '
   if (sgnum.eq. 63) sgname='Amma    '
   if (sgnum.eq. 63) sgname='Bbmm    '
   if (sgnum.eq. 63) sgname='Bmmb    '
   if (sgnum.eq. 63) sgname='Ccmm    '
   if (sgnum.eq. 63) sgname='Amam    '
   if (sgnum.eq. 64) sgname='Cmca    '
   if (sgnum.eq. 64) sgname='Abma    '
   if (sgnum.eq. 64) sgname='Bbcm    '
   if (sgnum.eq. 64) sgname='Bmab    '
   if (sgnum.eq. 64) sgname='Ccmb    '
   if (sgnum.eq. 64) sgname='Acam    '
   if (sgnum.eq. 65) sgname='Cmmm    '
   if (sgnum.eq. 65) sgname='Ammm    '
   if (sgnum.eq. 65) sgname='Bmmm    '
   if (sgnum.eq. 66) sgname='Cccm    '
   if (sgnum.eq. 66) sgname='Amaa    '
   if (sgnum.eq. 66) sgname='Bbmb    '
   if (sgnum.eq. 67) sgname='Cmma    '
   if (sgnum.eq. 67) sgname='Abmm    '
   if (sgnum.eq. 67) sgname='Bmcm    '
   if (sgnum.eq. 67) sgname='Bmam    '
   if (sgnum.eq. 67) sgname='Cmmb    '
   if (sgnum.eq. 67) sgname='Acmm    '
   if (sgnum.eq. 68) sgname='Ccca    '
   if (sgnum.eq. 68) sgname='Abaa    '
   if (sgnum.eq. 68) sgname='Bbcb    '
   if (sgnum.eq. 68) sgname='Bbab    '
   if (sgnum.eq. 68) sgname='Cccb    '
   if (sgnum.eq. 68) sgname='Acaa    '
   if (sgnum.eq. 69) sgname='Fmmm    '
   if (sgnum.eq. 70) sgname='Fddd    '
   if (sgnum.eq. 71) sgname='Immm    '
   if (sgnum.eq. 72) sgname='Ibam    '
   if (sgnum.eq. 72) sgname='Imcb    '
   if (sgnum.eq. 72) sgname='Icma    '
   if (sgnum.eq. 73) sgname='Ibca    '
   if (sgnum.eq. 73) sgname='Icab    '
   if (sgnum.eq. 74) sgname='Imma    '
   if (sgnum.eq. 74) sgname='Ibmm    '
   if (sgnum.eq. 74) sgname='Imcm    '
   if (sgnum.eq. 74) sgname='Imam    '
   if (sgnum.eq. 74) sgname='Immb    '
   if (sgnum.eq. 74) sgname='Icmm    '
   if (sgnum.eq. 75) sgname='P4      '
   if (sgnum.eq. 76) sgname='P41     '
   if (sgnum.eq. 77) sgname='P42     '
   if (sgnum.eq. 78) sgname='P43     '
   if (sgnum.eq. 79) sgname='I4      '
   if (sgnum.eq. 80) sgname='I41     '
   if (sgnum.eq. 81) sgname='P-4     '
   if (sgnum.eq. 82) sgname='I-4     '
   if (sgnum.eq. 83) sgname='P4/m    '
   if (sgnum.eq. 84) sgname='P42/m   '
   if (sgnum.eq. 85) sgname='P4/n    '
   if (sgnum.eq. 86) sgname='P42/n   '
   if (sgnum.eq. 87) sgname='I4/m    '
   if (sgnum.eq. 88) sgname='I41/a   '
   if (sgnum.eq. 89) sgname='P422    '
   if (sgnum.eq. 90) sgname='P4212   '
   if (sgnum.eq. 91) sgname='P4122   '
   if (sgnum.eq. 92) sgname='P41212  '
   if (sgnum.eq. 93) sgname='P4222   '
   if (sgnum.eq. 94) sgname='P42212  '
   if (sgnum.eq. 95) sgname='P4322   '
   if (sgnum.eq. 96) sgname='P43212  '
   if (sgnum.eq. 97) sgname='I422    '
   if (sgnum.eq. 98) sgname='I4122   '
   if (sgnum.eq. 99) sgname='P4mm    '
   if (sgnum.eq.100) sgname='P4bm    '
   if (sgnum.eq.101) sgname='P42cm   '
   if (sgnum.eq.102) sgname='P42nm   '
   if (sgnum.eq.103) sgname='P4cc    '
   if (sgnum.eq.104) sgname='P4nc    '
   if (sgnum.eq.105) sgname='P42mc   '
   if (sgnum.eq.106) sgname='P42bc   '
   if (sgnum.eq.107) sgname='I4mm    '
   if (sgnum.eq.108) sgname='I4cm    '
   if (sgnum.eq.109) sgname='I41md   '
   if (sgnum.eq.110) sgname='I41cd   '
   if (sgnum.eq.111) sgname='P-42m   '
   if (sgnum.eq.112) sgname='P-42c   '
   if (sgnum.eq.113) sgname='P-421m  '
   if (sgnum.eq.114) sgname='P-421c  '
   if (sgnum.eq.115) sgname='P-4m2   '
   if (sgnum.eq.116) sgname='P-4c2   '
   if (sgnum.eq.117) sgname='P-4b2   '
   if (sgnum.eq.118) sgname='P-4n2   '
   if (sgnum.eq.119) sgname='I-4m2   '
   if (sgnum.eq.120) sgname='I-4c2   '
   if (sgnum.eq.121) sgname='I-42m   '
   if (sgnum.eq.122) sgname='I-42d   '
   if (sgnum.eq.123) sgname='P4/mmm  '
   if (sgnum.eq.124) sgname='P4/mcc  '
   if (sgnum.eq.125) sgname='P4/nbm  '
   if (sgnum.eq.126) sgname='P4/nnc  '
   if (sgnum.eq.127) sgname='P4/mbm  '
   if (sgnum.eq.128) sgname='P4/mnc  '
   if (sgnum.eq.129) sgname='P4/nmm  '
   if (sgnum.eq.130) sgname='P4/ncc  '
   if (sgnum.eq.131) sgname='P42/mmc '
   if (sgnum.eq.132) sgname='P42/mcm '
   if (sgnum.eq.133) sgname='P42/nbc '
   if (sgnum.eq.134) sgname='P42/nnm '
   if (sgnum.eq.135) sgname='P42/mbc '
   if (sgnum.eq.136) sgname='P42/mnm '
   if (sgnum.eq.137) sgname='P42/nmc '
   if (sgnum.eq.138) sgname='P42/ncm '
   if (sgnum.eq.139) sgname='I4/mmm  '
   if (sgnum.eq.140) sgname='I4/mcm  '
   if (sgnum.eq.141) sgname='I41/amd '
   if (sgnum.eq.142) sgname='I41/acd '
   if (sgnum.eq.143) sgname='P3      '
   if (sgnum.eq.144) sgname='P31     '
   if (sgnum.eq.145) sgname='P32     '
   if (sgnum.eq.146) sgname='R3      '
   if (sgnum.eq.147) sgname='P-3     '
   if (sgnum.eq.148) sgname='R-3     '
   if (sgnum.eq.149) sgname='P312    '
   if (sgnum.eq.150) sgname='P321    '
   if (sgnum.eq.151) sgname='P3112   '
   if (sgnum.eq.152) sgname='P3121   '
   if (sgnum.eq.153) sgname='P3212   '
   if (sgnum.eq.154) sgname='P3221   '
   if (sgnum.eq.155) sgname='R32     '
   if (sgnum.eq.156) sgname='P3m1    '
   if (sgnum.eq.157) sgname='P31m    '
   if (sgnum.eq.158) sgname='P3c1    '
   if (sgnum.eq.159) sgname='P31c    '
   if (sgnum.eq.160) sgname='R3m     '
   if (sgnum.eq.161) sgname='R3c     '
   if (sgnum.eq.162) sgname='P-31m   '
   if (sgnum.eq.163) sgname='P-31c   '
   if (sgnum.eq.164) sgname='P-3m1   '
   if (sgnum.eq.165) sgname='P-3c1   '
   if (sgnum.eq.166) sgname='R-3m    '
   if (sgnum.eq.167) sgname='R-3c    '
   if (sgnum.eq.168) sgname='P6      '
   if (sgnum.eq.169) sgname='P61     '
   if (sgnum.eq.170) sgname='P65     '
   if (sgnum.eq.171) sgname='P62     '
   if (sgnum.eq.172) sgname='P64     '
   if (sgnum.eq.173) sgname='P63     '
   if (sgnum.eq.174) sgname='P-6     '
   if (sgnum.eq.175) sgname='P6/m    '
   if (sgnum.eq.176) sgname='P63/m   '
   if (sgnum.eq.177) sgname='P622    '
   if (sgnum.eq.178) sgname='P6122   '
   if (sgnum.eq.179) sgname='P6522   '
   if (sgnum.eq.180) sgname='P6222   '
   if (sgnum.eq.181) sgname='P6422   '
   if (sgnum.eq.182) sgname='P6322   '
   if (sgnum.eq.183) sgname='P6mm    '
   if (sgnum.eq.184) sgname='P6cc    '
   if (sgnum.eq.185) sgname='P63cm   '
   if (sgnum.eq.186) sgname='P63mc   '
   if (sgnum.eq.187) sgname='P-6m2   '
   if (sgnum.eq.188) sgname='P-6c2   '
   if (sgnum.eq.189) sgname='P-62m   '
   if (sgnum.eq.190) sgname='P-62c   '
   if (sgnum.eq.191) sgname='P6/mmm  '
   if (sgnum.eq.192) sgname='P6/mcc  '
   if (sgnum.eq.193) sgname='P63/mcm '
   if (sgnum.eq.194) sgname='P63/mmc '
   if (sgnum.eq.195) sgname='P23     '
   if (sgnum.eq.196) sgname='F23     '
   if (sgnum.eq.197) sgname='I23     '
   if (sgnum.eq.198) sgname='P213    '
   if (sgnum.eq.199) sgname='I213    '
   if (sgnum.eq.200) sgname='Pm-3    '
   if (sgnum.eq.201) sgname='Pn-3    '
   if (sgnum.eq.202) sgname='Fm-3    '
   if (sgnum.eq.203) sgname='Fd-3    '
   if (sgnum.eq.204) sgname='Im-3    '
   if (sgnum.eq.205) sgname='Pa-3    '
   if (sgnum.eq.206) sgname='Ia-3    '
   if (sgnum.eq.207) sgname='P432    '
   if (sgnum.eq.208) sgname='P4232   '
   if (sgnum.eq.209) sgname='F432    '
   if (sgnum.eq.210) sgname='F4132   '
   if (sgnum.eq.211) sgname='I432    '
   if (sgnum.eq.212) sgname='P4332   '
   if (sgnum.eq.213) sgname='P4132   '
   if (sgnum.eq.214) sgname='I4132   '
   if (sgnum.eq.215) sgname='P-43m   '
   if (sgnum.eq.216) sgname='F-43m   '
   if (sgnum.eq.217) sgname='I-43m   '
   if (sgnum.eq.218) sgname='P-43n   '
   if (sgnum.eq.219) sgname='F-43c   '
   if (sgnum.eq.220) sgname='I-43d   '
   if (sgnum.eq.221) sgname='Pm-3m   '
   if (sgnum.eq.222) sgname='Pn-3n   '
   if (sgnum.eq.223) sgname='Pm-3n   '
   if (sgnum.eq.224) sgname='Pn-3m   '
   if (sgnum.eq.225) sgname='Fm-3m   '
   if (sgnum.eq.226) sgname='Fm-3c   '
   if (sgnum.eq.227) sgname='Fd-3m   '
   if (sgnum.eq.228) sgname='Fd-3c   '
   if (sgnum.eq.229) sgname='Im-3m   '
   if (sgnum.eq.230) sgname='Ia-3d   '

   if (sgname.eq.'        ') then
      write(*,'(a,i4)') 'unknown space group number:',sgnum
      stop
   endif

 end subroutine getsgname

