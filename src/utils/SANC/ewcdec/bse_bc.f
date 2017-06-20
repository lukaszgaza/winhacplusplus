      SUBROUTINE BSE_BC_Wdec(Pg0bos,SWWbos,DRbos,SZZ_bd)   
*-
      IMPLICIT NONE
*-
      REAL*8 tHmu2,DREAL
      REAL*8 zm,rwm2,rzm2,rhm2
      REAL*8 ctw2,stw2,ctw4,ctw6
      REAL*8 rtz,rtw,rhw2,rhz2,rtw2,rhz,rhw,rth,rbw,rbw2,drtwbw
      REAL*8 rmf,Qf,cf,I3f
      REAL*8 Lnmuwm,Lnmuzm,Lnmuhm,Lnmutm,Lnmubm
      REAL*8 Pg0bos,SWWbos,DRbos,SZZ_bd
      REAL*8 b0fwwz,b0fzww,b0fwwh,b0fzzh,b0fww0
      REAL*8 b0pzww,b0pzhz
*-
      COMPLEX*16 cz2,wm2,zm2,hm2
      COMPLEX*16 b0f_Wdec,b0p_Wdec
*-
      COMMON/tHScal/tHmu2
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
      COMMON/BOS_m2/zm,cz2,wm2,zm2,hm2,rwm2,rzm2,rhm2
      COMMON/c_stw2/ctw2,stw2          
      COMMON/c_tw46/ctw4,ctw6             
      COMMON/RATIOS/rtz,rtw,rhw2,rhz2,rtw2,rhz,rhw,rth,rbw,rbw2,drtwbw
      COMMON/SCALES/Lnmuwm,Lnmuzm,Lnmuhm,Lnmutm,Lnmubm
*-
      b0fwwz = DREAL(b0f_Wdec(-rwm2,tHmu2,wm2,zm2))
      b0fwwh = DREAL(b0f_Wdec(-rwm2,tHmu2,wm2,hm2))
      b0fww0 = DREAL(b0f_Wdec(-rwm2,tHmu2,wm2,cz2))
*-
      b0fzww = DREAL(b0f_Wdec(-rzm2,tHmu2,wm2,wm2)) 
      b0fzzh = DREAL(b0f_Wdec(-rzm2,tHmu2,zm2,hm2))
      b0pzww = DREAL(b0p_Wdec(-rzm2,wm2,wm2))
      b0pzhz = DREAL(b0p_Wdec(-rzm2,hm2,zm2))
*-
*-
*---- SE BLOCK 
*-
      Pg0bos=2d0/3d0-3d0*Lnmuwm
*-
      SWWbos = + (1d0/12/ctw4+4d0/3/ctw2-17d0/3-4d0*ctw2)*b0fwwz
     &         + (1d0-1d0/3*rhw+1d0/12*rhw2)*b0fwwh
     &         - 4d0*stw2*b0fww0 
     &         + (3d0/2/rhw/ctw4+1d0/12/ctw4+2d0/3/ctw2-2d0)*Lnmuzm
     &         + (3d0/rhw-1d0/12/ctw2-1d0-1d0/12*rhw)*Lnmuwm
     &         + 1d0/12*rhw*(6d0+rhw)*Lnmuhm
     &         - (1d0/2/ctw4+1d0)/rhw-1d0/12/ctw4-3d0/4/ctw2
     &         - 1d0/9 - 1d0/12*rhw*(7d0+rhw)
*-
      DRbos=(1d0/12d0/ctw4+4d0/3d0/ctw2-17d0/3d0-4d0*ctw2)
     &      *(b0fwwz-ctw2*b0fzww) 
     &      +(1d0-1d0/3d0*rhw+1d0/12d0*rhw2)*b0fwwh
     &      -(1d0-1d0/3d0*rhz+1d0/12d0*rhz2)/ctw2*b0fzzh-4d0*stw2*b0fww0  
     &      -1d0/12d0*(-(1d0/ctw4+6d0/ctw2-24d0+rhw)*Lnmuzm
     &      -rhw2*stw2*Lnmuhm
     &      +(14d0+1d0/ctw2+16d0*ctw2-48d0*ctw4+rhw)*Lnmuwm
     &      +1d0/ctw4+19d0/3d0/ctw2-22d0/3d0+stw2*rhw2)
*-
      SZZ_bd=-(1d0/12d0-1d0/3d0*ctw2-3d0*ctw4)*b0fzww 
     &       -1d0/6d0*rhz*(1d0-1d0/2d0*rhz)*b0fzzh
     &       +(1d0/12d0+4d0/3d0*ctw2-17d0/3d0*ctw4
     &       -4d0*ctw6)*rzm2*b0pzww  
     &       +(1d0-1d0/3d0*rhz+1d0/12d0*rhz2)*rzm2*b0pzhz 
     &       +1d0/12d0*(1d0-rhz)*(Lnmuzm-rhz*Lnmuhm)
     &       -7d0/36d0+2d0/9d0*ctw2+1d0/6d0*rhz-1d0/12d0*rhz2
*-
      RETURN
      END
