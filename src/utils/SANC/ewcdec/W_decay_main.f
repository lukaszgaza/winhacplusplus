      PROGRAM  WDECAY
*--------------------
* Project SANC, an extraction from a manually written 
*               code "eeffLib" v0.20
* Authors: D.Yu. Bardin and L.V. Kalinovskaya
*-------------------
      IMPLICIT NONE!
*-------------------
      INTEGER*4 ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      INTEGER*4 VPOLTR,EWFFTR,GAMZTR,FERMTR,EWRCTR,EWWFFV,EMASTR
      INTEGER*4 TWOPOS,CHANNL,CHBLCK
      INTEGER*4 imu
*-
      REAL*8 tHmu2,tHla2,chfm,chfmp,wm,zm,hm
      REAL*8 omega,DSQRT
      REAL*8 DelVlim,DelV,DelSlim,DelS,DelH
      REAL*8 Born3,Virt3,Soft3,Hard3,TotWW,TotQED
      COMPLEX*16 FW_L,FW_R,FW_LD,FW_RD
*-
      COMMON/FLAGGS/ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      COMMON/TREATM/VPOLTR,EWFFTR,GAMZTR,FERMTR,EWRCTR,EWWFFV,EMASTR
      COMMON/HISREM/TWOPOS
      COMMON/CHANNL/CHANNL,CHBLCK
*-
*- Fixed flaggs
*-
      TWOPOS=1
      WEAK  =1
      FERMTR=1
*-
*- Channel flag setting
*-                   W decay
*- CHANNL=11--  electron--neutrino 
*- CHANNL=12--  muon    --neutrino
*- CHANNL=13--  tauon   --neutrino
*-
      IMOMS=1 ! =0 \alpha-scheme;  =1 GFermi-scheme
*-
      EWRCTR=2! =0 EWFF with QED contributions
*-            ! =1 EWFF with EW  contributions
*-            ! =2 EWFF with all contributions
*- 
      zm=91.1867d0
      hm=120d0
      wm=80.4514958d0 
*-      
      tHla2=1d0     !  mgamma^2
      omega=1d-10   !  ksoft
      PRINT *
      PRINT *,'mgamma^2=',tHla2,';  ksoft=',omega
*-
       DO 150 CHANNL=11,13
       IF    (CHANNL.EQ.11) THEN
        chfm =.51099907d-3
        chfmp=1d-10
        PRINT *
        PRINT *,'e-nu channel'
       ELSEIF(CHANNL.EQ.12) THEN
        chfm =.105658389d0
        chfmp=1d-10
        PRINT *
        PRINT *,'mu-nu channel'
       ELSEIF(CHANNL.EQ.13) THEN
        chfm =1.77705d0
        chfmp=1d-10
        PRINT *
        PRINT *,'tau-nu channel'
       ENDIF
*-
       DO imu=2,2
        tHmu2=(wm*(10d0)**(imu-2))**2
        PRINT *,'ch,masses=',CHANNL,chfm,chfmp
        DO EWRCTR=0,2,2
         PRINT 1,IMOMS,EWRCTR
 1       FORMAT(1x,'IMOMS,EWRCTR',2(3x,I2))
         CALL INEETT_Wdec(tHmu2,tHla2,chfm,chfmp,wm,zm,hm)
         CALL W_dec(omega,FW_L,FW_R,FW_LD,FW_RD,
     &              DelVlim,DelV,DelSlim,DelS,DelH,
     &              Born3,Virt3,Soft3,Hard3,TotWW,TotQED)
         PRINT 2,DSQRT(tHmu2),FW_L,FW_RD
 2       FORMAT(1x,'tHmu,FF_L,LD  =',F7.2,1x,2(F15.12,1x,F15.12,1x))
         PRINT 3,DelVlim
 3       FORMAT(1x,'DelVlim=',F13.10)
         PRINT 4,DelV
 4       FORMAT(1x,'DelV   =',F13.10)
         PRINT 5,DelSlim
 5       FORMAT(1x,'DelSlim=',F13.10)
         PRINT 6,DelS
 6       FORMAT(1x,'DelS   =',F13.10)
         PRINT 7,DelH
 7       FORMAT(1x,'DelH   =',F13.10)
         PRINT 8,DSQRT(tHmu2),Born3,Virt3,Soft3,Hard3,TotWW,TotQED
 8       FORMAT(1x,'W(GeV),BVSHTL=',F7.2,1x,3(F13.10,1x,F13.10,1x))
        ENDDO
       ENDDO
 150  CONTINUE
*-
      STOP
      END
