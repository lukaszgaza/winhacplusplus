      SUBROUTINE DELTR_Wdec(Deltar,Dr_bos,Dr_fer)
*-
      IMPLICIT NONE
*-
      INTEGER*4 ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      REAL*8 PI,PI2,D2,D3,D5,CEILER
      REAL*8 CONHC,alphaI,GFermi,GammaZ,DAL5H
      REAL*8 tHmu2,zm,rwm2,rzm2,rhm2
      REAL*8 ctw2,stw2,ctw4,ctw6,Lnctw
      REAL*8 rtz,rtw,rhw2,rhz2,rtw2,rhz,rhw,rth,rbw,rbw2,drtwbw
      REAL*8 Lnmuwm,Lnmuzm,Lnmuhm,Lnmutm,Lnmubm
      REAL*8 Pg0fer,DRfer,DRbos,DRwfer,DRwbos,Deltar,Dr_bos,Dr_fer,DREAL
*-
      COMPLEX*16 cz2,wm2,zm2,hm2
      COMPLEX*16 b0f_Wdec,b0fwwz,b0fzww,b0fwwh,b0fzzh,b0fww0
      COMPLEX*16 SWW_f_Wdec,SZZ_f_Wdec,Pgg_f_Wdec
*-
      COMMON/FLAGGS/ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      COMMON/MATHCO/PI,PI2,D2,D3,D5,CEILER
      COMMON/PHYSCO/CONHC,alphaI,GFermi,GammaZ,DAL5H
      COMMON/BOS_m2/zm,cz2,wm2,zm2,hm2,rwm2,rzm2,rhm2
      COMMON/c_stw2/ctw2,stw2          
      COMMON/c_tw46/ctw4,ctw6
      COMMON/RATIOS/rtz,rtw,rhw2,rhz2,rtw2,rhz,rhw,rth,rbw,rbw2,drtwbw
      COMMON/tHScal/tHmu2
      COMMON/SCALES/Lnmuwm,Lnmuzm,Lnmuhm,Lnmutm,Lnmubm
*-
*- Calculation of \Delta r in one-to-one with (6.390)-(6.395) of Thebook
*-
      b0fwwz=DREAL(b0f_Wdec(-rwm2,tHmu2,wm2,zm2))
      b0fwwh=DREAL(b0f_Wdec(-rwm2,tHmu2,wm2,hm2))
      b0fzww=DREAL(b0f_Wdec(-rzm2,tHmu2,wm2,wm2))
      b0fzzh=DREAL(b0f_Wdec(-rzm2,tHmu2,zm2,hm2))
      b0fww0=DREAL(b0f_Wdec(-rwm2,tHmu2,wm2,cz2))
*-
      Pg0fer=DREAL(Pgg_f_Wdec(0d0))
      Drfer =DREAL(SWW_f_Wdec(rwm2)-SZZ_f_Wdec(rzm2))/rwm2 
      Drwfer=DREAL(SWW_f_Wdec(0d0) -SWW_f_Wdec(rwm2))/rwm2
*- Formula (6.393)
*     DRbos =(1d0/12d0/ctw4+4d0/3d0/ctw2-17d0/3d0-4d0*ctw2)
*    &      *(b0fwwz-ctw2*b0fzww)
*    &      +(1d0-1d0/3d0*rhw+1d0/12d0*rhw2)*b0fwwh
*    &      -(1d0-1d0/3d0*rhz+1d0/12d0*rhz2)/ctw2*b0fzzh       
*    &      +1d0/12d0*stw2*rhw2*(Lnrhw-1d0)
*    &      -1d0/12d0*(1d0/ctw4+6d0/ctw2-24d0+rhw)*Lnctw
*    &      -1d0/12d0/ctw4-19d0/36d0/ctw2-133d0/18d0+8d0*ctw2
*- Representation coded in FF_lib, both are numerically equal (gauge invariance)
      DRbos =(1d0/12d0/ctw4+4d0/3d0/ctw2-17d0/3d0-4d0*ctw2)
     &      *(b0fwwz-ctw2*b0fzww)
     &      +(1d0-1d0/3d0*rhw+1d0/12d0*rhw2)*b0fwwh
     &      -(1d0-1d0/3d0*rhz+1d0/12d0*rhz2)/ctw2*b0fzzh-4d0*stw2*b0fww0
     &      -1d0/12d0*(-(1d0/ctw4+6d0/ctw2-24d0+rhw)*Lnmuzm
     &      -rhw2*stw2*Lnmuhm
     &      +(14d0+1d0/ctw2+16d0*ctw2-48d0*ctw4+rhw)*Lnmuwm
     &      +1d0/ctw4+19d0/3d0/ctw2-22d0/3d0+stw2*rhw2)
*-
      DRwbos=-(1d0/12d0/ctw4+4d0/3d0/ctw2-17d0/3d0-4d0*ctw2)*b0fwwz
     &      -(1d0-1d0/3d0*rhw+1d0/12d0*rhw2)*b0fwwh
     &      +3d0/4d0*rhw/(1d0-rhw)*(Lnmuhm-Lnmuwm)
     &      +1d0/4d0*rhw*(1d0-1d0/3d0*rhw)*Lnmuhm
     &      +(1d0/4d0-3d0/stw2)*(Lnmuwm-Lnmuzm)
     &      -(1d0/12d0/ctw4+17d0/12d0/ctw2)*Lnmuzm
     &      +(1d0/12d0/ctw2+11d0/2d0+1d0/12d0*rhw-4d0*stw2)*Lnmuwm 
     &      +1d0/12d0/ctw4+11d0/8d0/ctw2+139d0/36d0-177d0/24d0*ctw2
     &      +5d0/8d0*ctw4-1d0/12d0*rhw*(7d0/2d0-rhw)
*-
      Lnctw =Lnmuwm-Lnmuzm
      Deltar=1d0/alphaI/4d0/PI/stw2*(stw2*(-2d0/3-Pg0fer)
     &      +ctw2/stw2*(DRfer+DRbos)+DRwfer+DRwbos
     &      +11d0/2-5d0/8*ctw2*(1d0+ctw2)+9d0/4d0*ctw2/stw2*Lnctw
     &      +(3d0-7d0*ctw2)*Lnmuwm)
      Dr_bos=1d0/alphaI/4d0/PI/stw2*(stw2*(-2d0/3)
     &      +ctw2/stw2*DRbos+DRwbos
     &      +11d0/2-5d0/8*ctw2*(1d0+ctw2)+9d0/4d0*ctw2/stw2*Lnctw
     &      +(3d0-7d0*ctw2)*Lnmuwm)
      Dr_fer=1d0/alphaI/4d0/PI/stw2*(stw2*(-Pg0fer)
     &      +ctw2/stw2*DRfer+DRwfer)
*-
      return
      end
