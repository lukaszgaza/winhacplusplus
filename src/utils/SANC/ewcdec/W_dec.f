      SUBROUTINE W_dec(omega,F_L,F_R,F_LD,F_RD,
     &                 DelVlim,DelV,DelSlim,DelS,DelH,
     &                 Born3,Virt3,Soft3,Hard3,TotWW,TotQED)
*-----------------------------------------------------------
      IMPLICIT NONE!
*-------------------
      INTEGER*4 ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      INTEGER*4 VPOLTR,EWFFTR,GAMZTR,FERMTR,EWRCTR,EWWFFV,EMASTR
      INTEGER*4 CHANNL,CHBLCK,IHAMP,ZORO
*-
      REAL*8 omega,tHmu2,tHla2,reW,betae,lnreW
      REAL*8 sqrl,ruW,rdW,betaw,betau,betad,sqrle,betaew,betaeu,betaed
      REAL*8 PI,PI2,D2,D3,D5,CEILER
      REAL*8 CONHC,alphaI,GFermi,GammaZ,DAL5H
      REAL*8 ctw2,stw2,ctw4,ctw6
      REAL*8 wm,zm,tm,bm,rwm2,rzm2,rhm2,rtm2,rbm2,um,dm,rum2,rdm2,betaud
      REAL*8 qe,qt,qb,qe2,qt2,qb2,qeqt,anu,vnu,ae,ve,at,vt,ab,vb
      REAL*8 qtm,ae2,ve2,at2,vt2,ab2,vb2,vpae,vmae,vpat,vmat,vpab,vmab
      REAL*8 de,dt,Nc,coeff
      REAL*8 Lnmuwm,Lnmuzm,Lnmuhm,Lnmutm,Lnmubm,sdet3rw
      REAL*8 rtz,rtw,rhw2,rhz2,rtw2,rhz,rhw,rth,rbw,rbw2,drtwbw
      REAL*8 b0pcwt,b0ptwc,b0ptzt,b0pwwz,b0pwwh,b0ptht
      REAL*8 b0fwwh,b0fwwc,b0fwwz,b0fzzh,b0fzww
      REAL*8 Pg0fer,PZgzfer,SWWfer,DRfer,Pg0bos,SWWbos,DRbos,SZZ_bd
      REAL*8 BCTfer,BCTbos
      REAL*8 SWW_fd,SZZ_fd,SHH_fd,Tadfer,TotWST(0:1),TotWHA(0:1)
      REAL*8 Eb,Et,P,Nvu,sr2,FINT(12)
      REAL*8 DelSFT,ReJINTu,ReJINTd,ReC0IRF
      REAL*8 Delt_r,Dr_bos,Dr_fer
      REAL*8 sprmin,sprmax,xprmin,xprmax
      REAL*8 Born3,Virt3,Soft3,Hard3,TotWW,TotQED
      REAL*8 DSQRT,DLOG,DREAL,DDILOG_Wdec
      REAL*8 DelVlim,DelV,DelSlim,DelS,DelH
*-
      COMPLEX*16 c0_Wdec,ctwwcz,ctwtzc,ctwztw,ctwhtw
      COMPLEX*16 ctwt0b,ctwwb0,ctw0tw,JINTEG,JINTW_Wdec,C0IRF_Wdec,
     &           C0IR_Wdec
      COMPLEX*16 cz2,wm2,zm2,hm2,tm2,bm2,pm2
      COMPLEX*16 b0p_Wdec,b0f_Wdec,b0fwtc,b0ftwc,b0ftht,b0ftzt
      COMPLEX*16 bd1tbw,bd1ttz,bd1tth,bd2ttz,bdptbw,bdpttz
      COMPLEX*16 FFLQED,FFLEW,FFLD,F_L,F_R,F_LD,F_RD,HelAmp(4)
      COMPLEX*16 LQED,LCT,LP1,LP2,DCMPLX,DCONJG
*-
      COMMON/tHScal/tHmu2
      COMMON/tHScIR/tHla2
      COMMON/FLAGGS/ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      COMMON/TREATM/VPOLTR,EWFFTR,GAMZTR,FERMTR,EWRCTR,EWWFFV,EMASTR
      COMMON/CHANNL/CHANNL,CHBLCK
      COMMON/MATHCO/PI,PI2,D2,D3,D5,CEILER
      COMMON/PHYSCO/CONHC,alphaI,GFermi,GammaZ,DAL5H
      COMMON/SCALES/Lnmuwm,Lnmuzm,Lnmuhm,Lnmutm,Lnmubm
      COMMON/RATIOS/rtz,rtw,rhw2,rhz2,rtw2,rhz,rhw,rth,rbw,rbw2,drtwbw
      COMMON/c_stw2/ctw2,stw2  
      COMMON/c_tw46/ctw4,ctw6   
      COMMON/t_mass/tm,rtm2,tm2
      COMMON/b_mass/bm,rbm2,bm2
      COMMON/BOS_m2/zm,cz2,wm2,zm2,hm2,rwm2,rzm2,rhm2
      COMMON/COUPLC/qe,qt,qb,qe2,qt2,qb2,qeqt,anu,vnu,ae,ve,at,vt,ab,vb,
     &qtm,ae2,ve2,at2,vt2,ab2,vb2,vpae,vmae,vpat,vmat,vpab,vmab,de,dt,Nc
      COMMON/SUBBFS/bd1tbw,bd1ttz,bd1tth,bd2ttz,bdptbw,bdpttz
*-
      sr2=DSQRT(2d0)
      wm =DSQRT(rwm2)
      IF(WEAK.EQ.0) GOTO 111   
      CALL 
     &   BSE_FC_Wdec(Pg0fer,PZgzfer,SWWfer,DRfer,SWW_fd,SZZ_fd,SHH_fd,
     &               Tadfer)
      CALL BSE_BC_Wdec(Pg0bos,SWWbos,DRbos,SZZ_bd)

      BCTfer=1d0/2*(SWW_fd-stw2*Pg0fer+ctw2/stw2*DRfer)
      BCTbos=1d0/2*(-stw2*Pg0bos+ctw2/stw2*(DRbos-4d0*stw2*Lnmuwm))
*-                                          special contribution
      ctwt0b =C0IR_Wdec(tHla2,bm2,tm2,-rwm2)
      ctwwb0 =C0IR_Wdec(tHla2,bm2,wm2,-rtm2)
      ctw0tw =C0IR_Wdec(tHla2,tm2,wm2,-rbm2)
      ctwwcz=c0_Wdec(-rtm2,0d0,-rwm2,wm2,cz2,zm2)
      ctwtzc=c0_Wdec(-rtm2,0d0,-rwm2,tm2,zm2,cz2)
      ctwztw=c0_Wdec(-rtm2,0d0,-rwm2,zm2,tm2,wm2)
      ctwhtw=c0_Wdec(-rtm2,0d0,-rwm2,hm2,tm2,wm2)
*-
      b0ftwc=b0f_Wdec(-rtm2,tHmu2,wm2,cz2)
      b0ftzt=b0f_Wdec(-rtm2,tHmu2,zm2,tm2)
      b0ftht=b0f_Wdec(-rtm2,tHmu2,hm2,tm2)   
      b0fwtc=b0f_Wdec(-rwm2,tHmu2,tm2,cz2)
*-
      b0fwwh=DREAL(b0f_Wdec(-rwm2,tHmu2,wm2,hm2))
      b0fwwc=DREAL(b0f_Wdec(-rwm2,tHmu2,wm2,cz2))
      b0fwwz=DREAL(b0f_Wdec(-rwm2,tHmu2,wm2,zm2))
      b0fzzh=DREAL(b0f_Wdec(-rzm2,tHmu2,zm2,hm2))
      b0fzww=DREAL(b0f_Wdec(-rzm2,tHmu2,wm2,wm2))
*-
      b0ptzt=DREAL(b0p_Wdec(-rtm2,zm2,tm2))
      b0pwwz=DREAL(b0p_Wdec(-rwm2,wm2,zm2))
      b0pcwt=DREAL(b0p_Wdec(0d0,wm2,tm2))
      b0ptwc=DREAL(b0p_Wdec(-rtm2,wm2,cz2))
      b0pwwz=DREAL(b0p_Wdec(-rwm2,wm2,zm2))
      b0pwwh=DREAL(b0p_Wdec(-rwm2,wm2,hm2))
      b0ptht=DREAL(b0p_Wdec(-rtm2,hm2,tm2))
*-
      sdet3rw=rtw-1d0
*-
*-----------------------------------------------------------------------
*- W^+ --> u + d (dm-->0)
*-
      FFLQED = stw2*(DLOG(rwm2/tHla2)+qt2*DLOG(rtm2/tHla2)
     &                               +qb2*DLOG(rbm2/tHla2)
     &  -qt*qb*(2d0*(rtw+rbw-1d0)*rwm2*ctwt0b+3d0/2*Lnmutm+3d0/2*Lnmubm) 
     &  +qb   *(2d0*(1d0+rbw-rtw)*rwm2*ctwwb0+3d0/2*Lnmubm)
     &  -qt   *(2d0*(1d0+rtw-rbw)*rwm2*ctw0tw+3d0/2*Lnmutm))
*-
      FFLEW  = 
     & + stw2*qb2*(-(1d0-2d0/sdet3rw)*(b0fwtc-2d0)
     &                +2d0*rtw/sdet3rw*Lnmutm-4d0)
     & +1d0/2*(-vpab/ctw2*(1d0-2d0*ctw2-2d0*ctw4+2d0*ctw2*rtz 
     &         -(1d0-2d0*ctw2)/sdet3rw-2d0/sdet3rw**2)*rzm2*ctwwcz
     & +vpat*vpab/ctw2*(1d0+ctw2-rtz-1d0/sdet3rw+1d0/ctw2/sdet3rw**2)
     &                                                *rzm2*ctwtzc
     & -2d0*(vmat*ctw2*(1d0+rtw)+2d0*vpab+5d0*ctw2-1d0/4*rtw)
     &                                                *rzm2*ctwztw
     & +2d0*(1d0-1d0/4*rhw)                           *rtm2*ctwhtw
*-
     & +1d0/ctw2*(1d0/2*(vpat-vpab*vmab)*(1d0-2d0/sdet3rw)
     &                -vpab*vpat/ctw2/sdet3rw**2)*b0fwtc
     & +(vpab/ctw2*(1d0-1d0/sdet3rw-2d0/sdet3rw**2)+1d0/4-2d0/sdet3rw)
     &                                           *b0ftwc
     & +(vmab*vmat/ctw4*((1d0/2-ctw2)/sdet3rw+1d0/sdet3rw**2)
     &     -2d0*vmab/ctw2*(3d0/4+(1d0+ctw2)/sdet3rw)
     &     +19d0/8/ctw2-4d0+(1d0/2/ctw4+3d0/4/ctw2-1d0-4d0*ctw2)/sdet3rw
     &     +(1d0-2d0*ctw2)/ctw4/sdet3rw**2)      *b0ftzt
     & +(1d0-3d0/8*rhw+(1d0-1d0/4*rhw)/sdet3rw)  *b0ftht
     &  -(-1d0/12/ctw4+3d0/12/ctw2-5d0/3-4d0*ctw2
     &     +(3d0/4/ctw2-5d0-4d0*ctw2)/sdet3rw
     &           -2d0*vpab/ctw2*rtw/sdet3rw**2)  *b0fwwz
     & -(1d0-1d0/12*rhw*(1d0+rhw)+(1d0-1d0/4*rhw)/sdet3rw)
     &                                           *b0fwwh
     &  -1d0/4*vmat**2/ctw2*rtz*bd2ttz   !  term with difference
*-
     & -1d0/ctw4*(vmab*(vpab-2d0+4d0*ctw2)*(1d0/2+rtz)
     &      +5d0/4-3d0*ctw2+2d0*ctw4
     &      +(7d0/4-6d0*ctw2+4d0*ctw4)*rtz)  *rwm2*b0ptzt 
     &                              +4d0*stw2*rwm2*b0pwwz
     & -1d0/4*(1d0-rtw)*(2d0+rtw)*rwm2*(b0pcwt+b0ptwc) 
     & +1d0/3*(1d0/4/ctw4+4d0/ctw2-29d0)     *rwm2*b0pwwz
     & +(1d0-1d0/3*rhw+1d0/12*rhw2)          *rwm2*b0pwwh
     & -1d0/4*rtw*(1d0-4d0*rth)              *rhm2*b0ptht
*-
     & +(-1d0/12/ctw2+1d0/6+4d0*ctw2-1d0/12*rhw
     &           -(5d0/4-4d0*ctw2)/sdet3rw)*Lnmuwm
     &   -1d0/12*rhw*(5d0/2-rhw)           *Lnmuhm
     & -(1d0/2*vmab*(vpab/ctw2-2d0-(vmab-2d0)/ctw4/sdet3rw)
     &   -1d0/12/ctw4-17d0/24/ctw2+8d0/3-ctw2
     &                 -1d0/2/ctw4/sdet3rw)*Lnmuzm
     & +rtw*(5d0/4-4d0*(vmab+ctw2))/sdet3rw*Lnmutm
*-
     & +1d0/2*vmab*(3d0/2/ctw2+7d0-((vmab-2d0)/ctw4-16d0)/sdet3rw)
     & -1d0/12/ctw4-3d0/2/ctw2+301d0/36-1d0/2*ctw2
     & +1d0/2*rtw+7d0/24*rhw-1d0/12*rhw2-1d0/2/ctw4/sdet3rw)
     & +BCTfer+BCTbos
*-----------------------------------------------------------------------
      FFLD = 1d0/2*rtw/sdet3rw*(
     &  -vpab/ctw4*(2d0-5d0*ctw2+rtz-2d0/sdet3rw-6d0/sdet3rw**2)
     &                                           *rwm2*ctwwcz
     &  -1d0/ctw2*(vmab*(1d0+rtw)+3d0+4d0*ctw2-3d0/2*rtw)
     &                                           *rwm2*ctwztw
     &  +(1d0-3d0/2*rtw)                         *rhm2*ctwhtw
     &  +vpab/ctw4*(vmab*(2d0+3d0/ctw2/sdet3rw**2)
     &  -2d0+5d0*ctw2-rtz-2d0/sdet3rw-3d0*(1d0/ctw2-2d0)/sdet3rw**2) 
     &                                           *rwm2*ctwtzc
     &  +1d0/ctw2*(-3d0/4+3d0*ctw2-4d0*ctw4+(1d0/4+8d0*ctw2)/sdet3rw
     &                      +6d0*vpab/sdet3rw**2)     *b0fwwz
     &  +(1d0+1d0/4*rhw*rtw/sdet3rw)                  *b0fwwh
     &  +1d0/2/ctw2*(vpab*(1d0+4d0/sdet3rw
     &                      +6d0*(1d0/ctw2-2d0)/sdet3rw**2)
     &  -vmab*vpab*(1d0+6d0/ctw2/sdet3rw**2)+2d0*ctw2)*b0fwtc
     &  +1d0/ctw2*(ctw2+2d0*vpab*(1d0-1d0/sdet3rw-3d0/sdet3rw**2))
     &                                                *b0ftwc
     &  -1d0/2/ctw2*(vpab**2/ctw2*
     &           (1d0-2d0*ctw2-3d0/sdet3rw-6d0/sdet3rw**2)
     &                +4d0*vpab*(2d0*stw2-3d0/sdet3rw**2)
     &    -1d0+16d0*ctw2-8d0*ctw4+(1d0/2+16d0*ctw2)/sdet3rw)
     &                                                *b0ftzt
     &  -(rtw+1d0/2*rhw+1d0/4*rhw/sdet3rw)            *b0ftht
     &  -2d0*stw2*qb2*(b0fwtc+Lnmutm-2d0)
     &  -1d0/2*(8d0*ctw2-9d0*rtw/sdet3rw)*Lnmuwm
     &  -1d0/4*rhw                       *Lnmuhm
     &  -1d0/2/ctw2*(4d0*vpab+1d0/2+4d0*ctw2
     &             +vpab**2/ctw2*(1d0-3d0/sdet3rw))*Lnmuzm
     &  +1d0/2*(vpab**2/ctw2+8d0*vpab+8d0*ctw2 
     &          -2d0*rtw-9d0*rtw/sdet3rw)*Lnmutm
     &  +1d0/2*vmab*(vpab+1d0)/ctw4*(1d0-2d0*ctw2-3d0/sdet3rw)
     &          -vmab/ctw4*(1d0-4d0*ctw2+8d0*ctw4-3d0/sdet3rw)
     &  +1d0/2/ctw4-11d0/4/ctw2+9d0
     &      +1d0/4*rhw+3d0/2*rtw-3d0/2/ctw4/sdet3rw)     
     &  +1d0/2*(b0ftwc+Lnmuwm-1d0)
     &  +1d0/4/ctw4*(vmat*(vmab+2d0*ctw2)+1d0)*(b0ftzt+Lnmuzm-1d0) 
*-----------------------------------------------------------------------
*-
 111  coeff=1d0/alphaI/4d0/PI/stw2
*-
*---- Final expressions for the scalar form factors of the W decay
*-
      DO ZORO=0,1
       IF(WEAK.EQ.0) THEN
        F_L = DCMPLX(ZORO,0d0)
        F_LD= DCMPLX(0d0,0d0)
       ELSEIF(WEAK.EQ.1.AND.EWRCTR.EQ.0) THEN
        F_L = ZORO+coeff*FFLQED
        F_LD= (0d0,0d0)
       ELSEIF(WEAK.EQ.1.AND.EWRCTR.EQ.1) THEN
        F_L = ZORO+coeff*FFLEW
        F_LD=     +coeff*FFLD/rtm2
       ELSEIF(WEAK.EQ.1.AND.EWRCTR.EQ.2) THEN
        IF(IMOMS.EQ.0) THEN
         F_L = ZORO+coeff*(FFLQED+FFLEW)
         F_LD=     +coeff*FFLD/rtm2
        ELSE
         CALL DELTR_Wdec(Delt_r,Dr_bos,Dr_fer)
         F_L = ZORO+coeff*(FFLQED+FFLEW)-Delt_r/2
         F_LD=     +coeff*FFLD/rtm2
        ENDIF
       ENDIF
       F_R =(0d0,0d0)
       F_RD=(0d0,0d0)
*-
*---- Redefinition for calculation of probabilities
*-
       um=tm
       dm=bm
       rum2=rtm2
       rdm2=rbm2
       IF(CHANNL.EQ.11.OR.CHANNL.EQ.12.OR.CHANNL.EQ.13) THEN
        F_RD=F_LD
        F_LD=(0d0,0d0)
        um=bm
        dm=tm
        rum2=rbm2
        rdm2=rtm2
       ENDIF
       betaud=DSQRT(1d0-2d0*(rdm2+rum2)/rwm2+(rdm2-rum2)**2/rwm2**2)
*- ST=standard
       TotWST(ZORO)=GFermi*wm**3/(6d0*sr2*PI)*betaud*Nc*DREAL(
     &  +(1d0-1d0/2*(rdm2+rum2)/rwm2-1d0/2*(rdm2-rum2)**2/rwm2**2)
     &                 *(F_L*DCONJG(F_L)+F_R*DCONJG(F_R))
     &   +6d0*dm*um/rwm2*F_L*DCONJG(F_R)
     &  +betaud**2*(rdm2*F_L*DCONJG(F_RD)
     &             +rum2*F_L*DCONJG(F_LD)
     &            +dm*um*F_R*DCONJG(F_RD+F_LD)
     &   -2d0*rdm2*rum2*F_LD*DCONJG(F_RD)
     &    +1d0/2*(rwm2-rdm2-rum2)*(rdm2*F_RD*DCONJG(F_RD)
     &                            +rum2*F_LD*DCONJG(F_LD)) ) )
*- HA=Helicity Amplitudes
       Et=(rwm2+rum2-rdm2)/2d0/wm
       Eb=(rwm2+rdm2-rum2)/2d0/wm
       P =DSQRT(Et**2-rum2)
       Nvu=DSQRT(rwm2-(dm+um)**2)
       HelAmp(1)=1d0/sr2/Nvu*(
     & +        ( P*dm + P*um + dm*Et - um*Eb )*F_L
     & +        ( P*dm + P*um - dm*Et + um*Eb )*F_R
     & + 2*P*dm*( P**2 + P*wm + Eb*Et - dm*um )*F_RD
     & + 2*P*um*( P**2 - P*wm + Eb*Et - dm*um )*F_LD)
       HelAmp(2)=1d0/sr2/Nvu*(
     & +        ( P*dm + P*um - dm*Et + um*Eb )*F_L
     & +        ( P*dm + P*um + dm*Et - um*Eb )*F_R
     & + 2*P*dm*( P**2 - P*wm + Eb*Et - dm*um )*F_RD
     & + 2*P*um*( P**2 + P*wm + Eb*Et - dm*um )*F_LD)
       HelAmp(3)=1d0/Nvu*(
     & -        ( P**2 - P*wm + Eb*Et - dm*um )*F_L
     & +        ( P**2 + P*wm + Eb*Et - dm*um )*F_R )
       HelAmp(4)=1d0/Nvu*(
     & +        ( P**2 + P*wm + Eb*Et - dm*um )*F_L
     & -        ( P**2 - P*wm + Eb*Et - dm*um )*F_R )
       TotWHA(ZORO)=0d0
       DO IHAMP=1,4
        TotWHA(ZORO)=TotWHA(ZORO)+HelAmp(IHAMP)*DCONJG(HelAmp(IHAMP))
       ENDDO
      ENDDO
*-
      Virt3=GFermi*wm/(6d0*sr2*PI)*betaud*Nc*(TotWHA(1)-TotWHA(0))
*-
      IF(CHANNL.EQ.11.OR.CHANNL.EQ.12.OR.CHANNL.EQ.13) THEN
*-
*- Lepton case:
*-
       reW  =rdm2/rwm2
       betae=(1d0-reW)/(1d0+reW)
       lnreW=DLOG(reW)
*-
       Born3=GFermi*wm**3/(12d0*sr2*PI)*(1d0-reW)**2*(2d0+reW)
*- QED-only
*      Virt3=Born3/alphaI/PI/2d0*(
*    &  + 1d0/2*lnreW**2/betae+lnreW/betae*DLOG(rwm2/tHla2)
*    &  - 3d0/2*DLOG(rdm2/tHmu2)+DLOG(rwm2/tHla2)+DLOG(rdm2/tHla2) )
*      print *,'Virt3=',Born3+Virt3
*-
       Soft3=Born3/alphaI/PI/2d0*(
     &  - DLOG(4*omega**2/tHla2)*(lnreW/betae+2d0)
     &  - lnreW/betae
     &  + 2d0-1d0/betae*(DDILOG_Wdec(1d0-reW)
     &  - DDILOG_Wdec(1d0-1/reW)) )
       DelVlim=1/alphaI/PI*((DLOG(wm/dm)-1d0)*DLOG(tHla2/rwm2)
     &                      +DLOG(wm/dm)**2+1d0/2*DLOG(wm/dm))
       DelV=Virt3/Born3-1d0
       DelSlim=1/alphaI/PI*((DLOG(wm/dm)-1d0)*DLOG(4d0*omega**2/tHla2)
     &                      -DLOG(wm/dm)**2+DLOG(wm/dm)+1d0-D2)
       DelS=Soft3/Born3
*-
       Hard3=wm**3*GFermi/(16d0*sr2*PI)/alphaI/PI*(
     &  + 53d0/9-31d0/2*reW+5d0*reW**2+83d0/18*reW**3
     &  + 2d0*lnreW*(1d0-17d0/6*reW-7d0/3*reW**2))
     &  + Born3/alphaI/PI*( -2d0*DLOG(1d0-reW)
     &  - DLOG(rwm2/4/omega**2)*(1d0+1d0/2/betae*lnreW)
     &  + 1d0/betae*(DDILOG_Wdec(reW)-D2) )
       DelH=Hard3/Born3
*-
       xprmin=2d0*wm*omega
       xprmax=rwm2-rdm2
       sprmax=rwm2-xprmin
*
*   1) 1d0/(S-SPR)
       FINT(1)=-DLOG(xprmin/xprmax)
*   2) SPR
       FINT(2)=1d0/2*sprmax**2-1d0/2*rdm2**2
*   3) 1d0 
       FINT(3)=sprmax-rdm2
*   4) 1d0/SPR
       FINT(4)=DLOG(sprmax/rdm2)
*   5) 1d0/SPR**2
       FINT(5)=-1d0/sprmax+1d0/rdm2
*   6) 1d0/SPR**3
       FINT(6)=-1d0/2/sprmax**2+1d0/2/rdm2**2
*   7) JBETA/(S-SPR)
       FINT(7)=-DLOG(rwm2/rdm2)*DLOG(xprmin/xprmax)
     &        -DDILOG_Wdec(xprmax/rwm2)+DDILOG_Wdec(xprmin/rwm2)
*   8) JBETA*SPR
       FINT(8)=1d0/2*sprmax**2*DLOG(sprmax/rdm2)
     &        +1d0/4*(-sprmax**2+rdm2**2)
*   9) JBETA=LOG(SPR/M22)
       FINT(9)=sprmax*DLOG(sprmax/rdm2)-sprmax+rdm2
*-
       Hard3=GFermi*wm/(36d0*sr2*PI)/alphaI/PI*(
     &  +2d0*rdm2**3                   *FINT(6)
     &  -2d0*rdm2**2*(3d0+reW)         *FINT(5) 
     &  +3d0/2*rdm2*(6d0+5d0*reW)      *FINT(4) 
     &  +(5d0+3d0/2*reW)/rwm2          *FINT(2) 
     &  +(7d0-33d0/2*reW-15d0/2*reW**2)*FINT(3)
     &  -3d0/2*(2d0+reW)/rwm2          *FINT(8)
     &  -3d0*(1d0-3d0/2*reW-reW**2)    *FINT(9))
     &  +Born3/alphaI/PI*(1d0/betae*FINT(7)-2d0*FINT(1))
*-
      TotQED=Born3*(1d0
     &      +1d0/alphaI/PI*(77d0/24-2d0*D2-3d0/4*DLOG(rwm2/tHmu2)))
*-
      ELSE
*-
*- ud(cs) case:
*-
       sqrl =DSQRT(rwm2**2-2d0*rwm2*(rtm2+rbm2)+(rtm2-rbm2)**2)
       betaw=sqrl/(rwm2-rtm2-rbm2)
       betau=sqrl/(rwm2+rtm2-rbm2)
       betad=sqrl/(rwm2-rtm2+rbm2)
       ruW=rtm2/rwm2
       rdW=rbm2/rwm2
*-
       Born3=GFermi*wm/(12d0*sr2*PI)*sqrl*Nc
     &      *((1d0-ruW+rdW)*(1d0+ruW-rdW)+1d0-ruW-rdW)
*-
       ReJINTu=DREAL(JINTW_Wdec(bm2,wm2,-rtm2))
       ReJINTd=DREAL(JINTW_Wdec(tm2,wm2,-rbm2))
       ReC0IRF=DREAL(C0IRF_Wdec(tm2,bm2,-rwm2))
       Soft3=Born3/alphaI/2d0/PI*(+2d0
* - assotiated with `-um^2', `qb' and `betad'
     & - qb*(qb/betad+(qt/betaw-1d0/betad)*DLOG(4d0*omega**2/tHla2))
     &                                    *ReJINTu
* - assotiated with `-dm^2', `qt' and `betau'
     & - qt*(qt/betau+(qb/betaw+1d0/betau)*DLOG(4d0*omega**2/tHla2))
     &                                    *ReJINTd
     & - (1d0+qt**2+qb**2)*DLOG(4d0*omega**2/tHla2)
     & + 2d0*qt*qb*((rwm2-rtm2-rbm2)*ReC0IRF-PI2/betaw)
     & + qt/betau*(DDILOG_Wdec(-1d0/2*((rwm2+rtm2-rbm2)**2*betau
     &                    +(rwm2+rtm2-rbm2)**2*betau**2)/rwm2/rtm2)
     &           - DDILOG_Wdec(2d0/(1d0+betau)*betau))
     & - qb/betad*(DDILOG_Wdec(-1d0/2*((rwm2-rtm2+rbm2)**2*betad
     &                    +(rwm2-rtm2+rbm2)**2*betad**2)/rwm2/rbm2)
     &           - DDILOG_Wdec(2d0/(1d0+betad)*betad)) )
*-
       xprmin=2d0*wm*omega
       sprmax=rwm2-xprmin
       sprmin=(tm+bm)**2
       xprmax=rwm2-sprmin  !  Penka, it is unused!!!
       sqrle =DSQRT(sprmax**2-2d0*sprmax*(rtm2+rbm2)+(rtm2-rbm2)**2)
       betaew=sqrle/(sprmax-rtm2-rbm2)
       betaeu=sqrle/(sprmax+rtm2-rbm2)
       betaed=sqrle/(sprmax-rtm2+rbm2)
*
*   1) sqrlpr/(s-spr)
       FINT(1)=-sqrle-1d0/2*(rwm2-rtm2-rbm2)
     &      *DLOG(((1d0+betaew)*(sprmax-rtm2-rbm2))**2/(4d0*rtm2*rbm2))
     & + sqrl*DLOG((sqrl+sqrle-xprmin)**2
     &   /2d0/xprmin/(rwm2-rtm2-rbm2+sqrle-xprmin))
     & - sqrl*DLOG((sqrl-rwm2+sprmin)/(sqrl+rwm2-sprmin))
*   2) sqrlpr*spr
       FINT(2)=sqrle**3/3d0+ sqrle*(rtm2+rbm2)*(sprmax-rtm2-rbm2)/2
     &        -rtm2*rbm2*(rtm2+rbm2)*DLOG((1d0+betaew)/(1d0-betaew))
*   3) sqrlpr 
       FINT(3)=sqrle*(sprmax-rtm2-rbm2)/2d0
     &                     - rtm2*rbm2*DLOG((1d0+betaew)/(1d0-betaew))
*   4) sqrlpr/spr
       FINT(4)=sqrle-1d0/2*(rtm2+rbm2)*DLOG((1d0+betaew)/(1d0-betaew))
     &              -1d0/2*(rtm2-rbm2)*DLOG((1d0+betaeu)/(1d0-betaeu))
     &              +1d0/2*(rtm2-rbm2)*DLOG((1d0+betaed)/(1d0-betaed))
*   5) sqrlpr/spr**2
       FINT(5)=-sqrle/sprmax  +1d0/2*DLOG((1d0+betaew)/(1d0-betaew))
     &    +1d0/2*(rtm2+rbm2)/(rtm2-rbm2)*DLOG((1d0+betaeu)/(1d0-betaeu))
     &    -1d0/2*(rtm2+rbm2)/(rtm2-rbm2)*DLOG((1d0+betaed)/(1d0-betaed))
*   6) sqrlpr/spr**3
       FINT(6)=sqrle/sprmax/2*(rtm2+rbm2)/(rtm2-rbm2)**2
     &        -sqrle/sprmax**2/2
     &       +rtm2*rbm2/(rtm2-rbm2)**3*DLOG((1d0+betaeu)/(1d0-betaeu))
     &       -rtm2*rbm2/(rtm2-rbm2)**3*DLOG((1d0+betaed)/(1d0-betaed))
*   7) JBETA1/(s-spr) - assotiated with `-dm^2', `qt' and `betau'
       FINT(7)=DLOG((1d0+betau)/(1d0-betau))*(
     &    -DLOG((sqrl+sqrle-xprmin)**2
     &                    /2d0/xprmin/(rwm2-rtm2-rbm2+sqrle-xprmin)) 
     &    +DLOG((sqrl-rwm2+sprmin)/(sqrl+rwm2-sprmin)) )
     &    +DDILOG_Wdec(2d0*xprmin*(rwm2-rtm2-rbm2+sqrle-xprmin)
     &        /(sqrl+sqrle-xprmin)/(sqrl+rwm2-rtm2+rbm2))
     &    -DDILOG_Wdec((sqrl+rwm2-sprmin)/(sqrl+rwm2-rtm2+rbm2))
     &    -DDILOG_Wdec(2d0*xprmin*(rwm2-rtm2-rbm2+sqrle-xprmin)
     &        /(sqrl+sqrle-xprmin)/(sqrl+rwm2-rtm2-rbm2))
     &    +DDILOG_Wdec((sqrl+rwm2-sprmin)/(sqrl+rwm2-rtm2-rbm2))
     &    -DDILOG_Wdec(2d0*xprmin*(rwm2-rtm2-rbm2+sqrle-xprmin)
     &        /(sqrl+sqrle-xprmin)/(sqrl+rwm2+rtm2-rbm2))
     &    +DDILOG_Wdec((sqrl+rwm2-sprmin)/(sqrl+rwm2+rtm2-rbm2))
     &    +DDILOG_Wdec((sqrl+sqrle-xprmin)/(sqrl-rwm2+rtm2-rbm2))
     &    -DDILOG_Wdec((sqrl-rwm2+sprmin)/(sqrl-rwm2+rtm2-rbm2))
     &    -DDILOG_Wdec((sqrl+sqrle-xprmin)/(sqrl-rwm2+rtm2+rbm2))
     &    +DDILOG_Wdec((sqrl-rwm2+sprmin)/(sqrl-rwm2+rtm2+rbm2))
     &    -DDILOG_Wdec((sqrl+sqrle-xprmin)/(sqrl-rwm2-rtm2+rbm2))
     &    +DDILOG_Wdec((sqrl-rwm2+sprmin)/(sqrl-rwm2-rtm2+rbm2))
     &    +DDILOG_Wdec((sqrle-sprmax+rtm2-rbm2)
     &                /(sqrle+sprmax+rtm2-rbm2))
     &    +DDILOG_Wdec((sqrle+sprmax-rtm2+rbm2)
     &                /(sqrle-sprmax-rtm2+rbm2))
     &    -2d0*DDILOG_Wdec(-DSQRT(rbm2/rtm2))
*   8) JBETA1*spr
       FINT(8)=1d0/4*sqrle*(sprmax+rtm2+5d0*rbm2)
     &    -1d0/4*(rtm2-rbm2)*(rtm2+rbm2)*DLOG((1d0+betaew)/(1d0-betaew))
     &    +1d0/4*(-2d0*sprmax**2+(rtm2+rbm2)**2+2d0*rtm2*rbm2)
     &                                  *DLOG((1d0+betaeu)/(1d0-betaeu))
     &    +1d0/4*((rtm2+rbm2)**2+2d0*rtm2*rbm2)
     &                                  *DLOG((1d0+betaed)/(1d0-betaed))
*  9) JBETA1
       FINT(9)=sqrle-1d0/2*(rtm2-rbm2)*DLOG((1d0+betaew)/(1d0-betaew))
     &    +(-sprmax+1d0/2*(rtm2+rbm2))*DLOG((1d0+betaeu)/(1d0-betaeu))
     &             +1d0/2*(rtm2+rbm2) *DLOG((1d0+betaed)/(1d0-betaed))
*  10) JBETA2/(s-spr) - assotiated with `-um^2', `qb' and `betad'
       FINT(10)=DLOG((1d0+betad)/(1d0-betad))*(
     &    -DLOG((sqrl+sqrle-xprmin)**2
     &                    /2d0/xprmin/(rwm2-rtm2-rbm2+sqrle-xprmin)) 
     &    +DLOG((sqrl-rwm2+sprmin)/(sqrl+rwm2-sprmin)) )
     &    +DDILOG_Wdec(2d0*xprmin*(rwm2-rtm2-rbm2+sqrle-xprmin)
     &        /(sqrl+sqrle-xprmin)/(sqrl+rwm2+rtm2-rbm2))
     &    -DDILOG_Wdec((sqrl+rwm2-sprmin)/(sqrl+rwm2+rtm2-rbm2))
     &    -DDILOG_Wdec(2d0*xprmin*(rwm2-rtm2-rbm2+sqrle-xprmin)
     &        /(sqrl+sqrle-xprmin)/(sqrl+rwm2-rtm2-rbm2))
     &    +DDILOG_Wdec((sqrl+rwm2-sprmin)/(sqrl+rwm2-rtm2-rbm2))
     &    -DDILOG_Wdec(2d0*xprmin*(rwm2-rtm2-rbm2+sqrle-xprmin)
     &        /(sqrl+sqrle-xprmin)/(sqrl+rwm2-rtm2+rbm2))
     &    +DDILOG_Wdec((sqrl+rwm2-sprmin)/(sqrl+rwm2-rtm2+rbm2))
     &    +DDILOG_Wdec((sqrl+sqrle-xprmin)/(sqrl-rwm2-rtm2+rbm2))
     &    -DDILOG_Wdec((sqrl-rwm2+sprmin)/(sqrl-rwm2-rtm2+rbm2))
     &    -DDILOG_Wdec((sqrl+sqrle-xprmin)/(sqrl-rwm2+rtm2+rbm2))
     &    +DDILOG_Wdec((sqrl-rwm2+sprmin)/(sqrl-rwm2+rtm2+rbm2))
     &    -DDILOG_Wdec((sqrl+sqrle-xprmin)/(sqrl-rwm2+rtm2-rbm2))
     &    +DDILOG_Wdec((sqrl-rwm2+sprmin)/(sqrl-rwm2+rtm2-rbm2))
     &  +DDILOG_Wdec((sqrle-sprmax-rtm2+rbm2)/(sqrle+sprmax-rtm2+rbm2))
     &  +DDILOG_Wdec((sqrle+sprmax+rtm2-rbm2)/(sqrle-sprmax+rtm2-rbm2))
     &    -2d0*DDILOG_Wdec(-DSQRT(rtm2/rbm2))
*  11) JBETA2*spr
       FINT(11)=1d0/4*sqrle*(sprmax+rbm2+5d0*rtm2)
     &    +1d0/4*(rtm2-rbm2)*(rtm2+rbm2)*DLOG((1d0+betaew)/(1d0-betaew))
     &    +1d0/4*((rtm2+rbm2)**2+2d0*rtm2*rbm2)
     &                                  *DLOG((1d0+betaeu)/(1d0-betaeu))
     &    +1d0/4*(-2d0*sprmax**2+(rtm2+rbm2)**2+2d0*rtm2*rbm2)
     &                                  *DLOG((1d0+betaed)/(1d0-betaed))
*  12) JBETA2
       FINT(12)=sqrle+1d0/2*(rtm2-rbm2)*DLOG((1d0+betaew)/(1d0-betaew))
     &               +1d0/2*(rtm2+rbm2)*DLOG((1d0+betaeu)/(1d0-betaeu))
     &               -1d0/2*(rtm2+rbm2)*DLOG((1d0+betaed)/(1d0-betaed))
     &              -(sprmax-rtm2-rbm2)*DLOG((1d0+betaed)/(1d0-betaed))
*---------------------------------------------------------------------->
       Hard3=GFermi*wm/(12d0*sr2*PI)/alphaI/PI*(
     &  - 2d0/3*(ruW-rdW)**2*rwm2**2                          *FINT(6)
     &  + 2d0/3*((3d0*qt-1d0)*ruW-(3d0*qb+1d0)*rdW+(ruW-rdW)**2)*rwm2
     &                                                        *FINT(5)
     &  -(2d0/3+qt2+qb2-(1d0/6-2d0*qt)*ruW-(1d0/6+2d0*qb)*rdW)*FINT(4)
     &  +(5d0/3+2d0*qt*qb+1d0/2*(ruW+rdW))/rwm2               *FINT(3)
     &  + qb2*(1d0+1d0/2*(ruW+rdW))/rwm2        *FINT(11)
     &  + qb2*(1d0-3d0/2*(ruW+rdW)-(ruW-rdW)**2)*FINT(12)
     &  + qt2*(1d0+1d0/2*(ruW+rdW))/rwm2        *FINT(8)
     &  + qt2*(1d0-3d0/2*(ruW+rdW)-(ruW-rdW)**2)*FINT(9)
     &                                         )
     &  + Born3/alphaI/PI*(
     &  - (1d0+qt2+qb2)/sqrl         *FINT(1)
     &  + qb*(2d0*rbm2/sqrl-qb/betaw)*FINT(10)
     &  - qt*(2d0*rtm2/sqrl+qt/betaw)*FINT(7)
     &                    )
      ENDIF
*-
      TotWW =Virt3+Soft3+Hard3
      TotQED=Born3*(1d0
     &      +1d0/alphaI/PI*(77d0/24-2d0*D2-3d0/4*DLOG(rwm2/tHmu2)))
*-
      RETURN
      END
