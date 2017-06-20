      SUBROUTINE INEETT_Wdec(tHmu2A,tHla2A,tmA,bmA,wmA,zmA,hmA)
*---------------------------------------------------------
      IMPLICIT NONE!
*-------------------
      INTEGER*4 ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      INTEGER*4 VPOLTR,EWFFTR,GAMZTR,FERMTR,EWRCTR,EWWFFV,EMASTR
      INTEGER*4 TWOPOS,CHANNL,CHBLCK
      INTEGER*4 imf
      REAL*8 PI,PI2,D2,D3,D5,CEILER
      REAL*8 CONHC,alphaI,GFermi,GammaZ,DAL5H
      REAL*8 tHmu2,tHmu2A,tHla2,tHla2A
      REAL*8 qe,qt,qb,qe2,qt2,qb2,qeqt,anu,vnu,ae,ve,at,vt,ab,vb,i3t,qtc
      REAL*8 qtm,ae2,ve2,at2,vt2,ab2,vb2,vpae,vmae,vpat,vmat,vpab,vmab
      REAL*8 vpan,vman,vpau,vmau,vpad,vmad
      REAL*8 i3e,de,dt,Nc
      REAL*8 rmf1,rmf2,rmf3,rmf,Qf,cf,I3f
      REAL*8 wm,wmA,zm,zmA,hm,hmA,tm,tmA,rwm2,rzm2,rhm2,rtm2
      REAL*8 ctw2,stw2,ctw4,ctw6
      REAL*8 rtz,rtw,rhw2,rhz2,rtw2,rhz,rhw,rth,rbw,rbw2,drtwbw
      REAL*8 em,rem2,bm,bmA,rbm2,lnrbw,lnrtz,lnrth 
      REAL*8 Lnmuwm,Lnmuzm,Lnmuhm,Lnmutm,Lnmubm
      REAL*8 DLOG,DABS
      COMPLEX*16 cz2,wm2,zm2,hm2,tm2,bm2,em2,
     &           b0f_Wdec,b0p_Wdec,b0d1_Wdec,
     &           b0d2_Wdec,b0dp_Wdec,b0dw_Wdec
      COMPLEX*16 bd1tbw,bd1ttz,bd1tth,bd2ttz,bdptbw,bdpttz,DCMPLX,test
      LOGICAL*4 FIRST
      SAVE FIRST
      DATA FIRST/.FALSE./
*-
      DIMENSION rmf1(12),rmf2(12),rmf3(12)
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
*-        nue,e,numu,mu,nutau,tau,um,dm,cm,sm,tm,bm
*-
      COMMON/FLAGGS/ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      COMMON/TREATM/VPOLTR,EWFFTR,GAMZTR,FERMTR,EWRCTR,EWWFFV,EMASTR
      COMMON/HISREM/TWOPOS
      COMMON/CHANNL/CHANNL,CHBLCK
      COMMON/MATHCO/PI,PI2,D2,D3,D5,CEILER
      COMMON/PHYSCO/CONHC,alphaI,GFermi,GammaZ,DAL5H
      COMMON/tHScal/tHmu2
      COMMON/tHScIR/tHla2
      COMMON/t_mass/tm,rtm2,tm2
      COMMON/b_mass/bm,rbm2,bm2
      COMMON/e_mass/em,rem2,em2
*-------------------------->em2 or cz2 depending on EMASTR(1,0)
      COMMON/BOS_m2/zm,cz2,wm2,zm2,hm2,rwm2,rzm2,rhm2
      COMMON/c_stw2/ctw2,stw2
      COMMON/c_tw46/ctw4,ctw6
      COMMON/RATIOS/rtz,rtw,rhw2,rhz2,rtw2,rhz,rhw,rth,rbw,rbw2,drtwbw
      COMMON/COUPLC/qe,qt,qb,qe2,qt2,qb2,qeqt,anu,vnu,ae,ve,at,vt,ab,vb,
     &qtm,ae2,ve2,at2,vt2,ab2,vb2,vpae,vmae,vpat,vmat,vpab,vmab,de,dt,Nc
      COMMON/COUPCC/vpan,vman,vpau,vmau,vpad,vmad
      COMMON/I3CHNL/i3t,qtc
      COMMON/SCALES/Lnmuwm,Lnmuzm,Lnmuhm,Lnmutm,Lnmubm
      COMMON/SUBBFS/bd1tbw,bd1ttz,bd1tth,bd2ttz,bdptbw,bdpttz
*-
      DATA rmf1/1d-10,.51099907d-3,1d-10,.105658389d0,1d-10,1.77705d0,
     &         .062d0,.083d0,1.50d0,.215d0,173.8d0,4.70d0/
      DATA rmf2/1d-10,.51099907d-3,1d-10,.105658389d0,1d-10,1.77705d0,
     &        .04851d0,.04851d0,1.50d0,.150d0,173.8d0,4.70d0/
*-             to reproduce ALQEDI=  128.88667 (Eidelman-Jegerlehner)
      DATA rmf3/1d-10,.51099907d-3,1d-10,.105658389d0,1d-10,1.77705d0,
     &        .075d0,.075d0,1.50d0,.250d0,173.8d0,4.70d0/
*-             to reproduce dispersive treatment of \Pi_{Z\gamma}(0)
      DATA Qf/0d0,-1d0,0d0,-1d0,0d0,-1d0,
     &        .666666666666667d0,-.333333333333333d0,
     &        .666666666666667d0,-.333333333333333d0,
     &        .666666666666667d0,-.333333333333333d0/
      DATA cf/6*1d0,6*3d0/
      DATA I3f/.5d0,-.5d0,.5d0,-.5d0,.5d0,-.5d0,.5d0,-.5d0,
     &         .5d0,-.5d0,.5d0,-.5d0/
***   DATA alphaI/137.035989 5 D0/,CONHC/.389 379 66 D+9/ ! Original setting.
      DATA alphaI/137.035999 76D0/,CONHC/.389 379 66 D+9/ ! For comp with FA.
*-                                 conversion factor to pb's
*---------------------------------------- initialization --------------------
*-
      PI =4D0*DATAN(1D0)
      PI2=PI**2
      D2 =PI2/6D0
      D3 =1.2020569031596D0
      D5 =1.0369277551434D0
      CEILER=.577216D0      
*-
*-  Masses
*-
*- One had to be more careful with these three lines of Zeuthen !!!
*- Sometimes they might be needed, however, and their commenting 
*- is not a satisfactory solution !!!
*-
*     rmf1(11)=tmA
*     rmf2(11)=tmA
*     rmf3(11)=tmA
*-
      IF(FERMTR.EQ.1) THEN
       DO imf=1,12
        rmf(imf)=rmf1(imf)
       ENDDO
      ELSEIF(FERMTR.EQ.2) THEN
       DO imf=1,12
        rmf(imf)=rmf2(imf)
       ENDDO
      ELSEIF(FERMTR.EQ.3) THEN
       DO imf=1,12
        rmf(imf)=rmf3(imf)
       ENDDO
      ENDIF
*-
      em=rmf(2)
      tm=tmA
      bm=bmA
      wm=wmA
      zm=zmA
      hm=hmA
*-  
      rem2=em**2
      rtm2=tm**2
      rbm2=bm**2
      rwm2=wm**2
      rzm2=zm**2
      rhm2=hm**2
*-
      cz2 =DCMPLX(0D0, -1D-30)
      IF(EMASTR.EQ.0) THEN
       em2=cz2
      ELSE            
       em2=DCMPLX(rem2,-1D-30)
      ENDIF
      tm2 =DCMPLX(rtm2,-1D-30)                 
      bm2 =DCMPLX(rbm2,-1D-30)                 
      wm2 =DCMPLX(rwm2,-1D-30)
      zm2 =DCMPLX(rzm2,-1D-30)
      hm2 =DCMPLX(rhm2,-1D-30)                 
*-
*-  Mass ratios 
*-
      ctw2=rwm2/rzm2
      stw2=1d0-ctw2     
      ctw4=ctw2**2
      ctw6=ctw2**3
      rtz =rtm2/rzm2
      rtw =rtm2/rwm2
      rtw2=rtw**2
      rth =rtm2/rhm2
      rhz =rhm2/rzm2
      rhz2=rhz**2
      rhw =rhm2/rwm2
      rhw2=rhw**2
      rbw =rbm2/rwm2
      rbw2=rbw**2
      drtwbw=rtw-rbw
*-
      IF(TWOPOS.EQ.1) THEN
       GammaZ=2.499 776 D0 ! First  IPS
      ELSE
       GammaZ=2.499 538 D0 ! Second IPS
      ENDIF
*-
*-  For comparison with Hollik
*-
      IF(GAMZTR.EQ.0) GammaZ=0d0
*-
      IF(IMOMS.EQ.0) THEN 
       GFermi=pi/alphaI/sqrt(2d0)/stw2/ctw2/rzm2
      ELSE
       GFermi=1.16637D-5
      ENDIF      
*-
*-  Scales
*-
      tHmu2  = tHmu2A
      tHla2  = tHla2A
      IF(rbm2.EQ.0d0) Lnmubm=0d0
      IF(rbm2.NE.0d0) Lnmubm = DLOG(rbm2/tHmu2)
      Lnmutm = DLOG(rtm2/tHmu2)
      Lnmuwm = DLOG(rwm2/tHmu2) 
      Lnmuzm = DLOG(rzm2/tHmu2)
      Lnmuhm = DLOG(rhm2/tHmu2)
*-
*-  Expansions obtained with 
*-
      bd1ttz=DCMPLX(0d0,0d0)
      bd2ttz=DCMPLX(0d0,0d0)
      bdpttz=DCMPLX(0d0,0d0)
      bd1tth=DCMPLX(0d0,0d0)
      bd1tbw=DCMPLX(0d0,0d0)
      bdptbw=DCMPLX(0d0,0d0)
      IF(CHANNL.EQ.9.AND.CHBLCK.EQ.0) GOTO 111
      IF(rtz.LE.1d-3) THEN
        lnrtz=DLOG(rtz)
        bd1ttz=lnrtz+1d0/2+rtz*(2d0*lnrtz+5d0/3)
     &                 +rtz**2*(5d0*lnrtz+59d0/12)
     &                 +rtz**3*(14d0*lnrtz+449d0/30)
     &                 +rtz**4*(42d0*lnrtz+1417d0/30)
     &                 +rtz**5*(132d0*lnrtz+16127d0/105)
        bd2ttz=2d0*lnrtz+5d0/3+rtz*(5d0*lnrtz+59d0/12)
     &                 +rtz**2*(14d0*lnrtz+449d0/30)
     &                 +rtz**3*(42d0*lnrtz+1417d0/30)
     &                 +rtz**4*(132d0*lnrtz+16127d0/105)
     &                 +rtz**5*(429d0*lnrtz+429697d0/840)
        bdpttz=-lnrtz-17d0/6-rtz*(7d0*lnrtz+133d0/12)
     &                 -rtz**2*(31d0*lnrtz+2657d0/60)
     &                 -rtz**3*(126d0*lnrtz+5231d0/30)
     &                 -rtz**4*(498d0*lnrtz+142991d0/210)
     &                 -rtz**5*(1947d0*lnrtz+2225231d0/840)
**     print *,'rtz=',rtz
**     print 1,bd1ttz
**     print 4,lnrtz+1d0/2
**     print 4,rtz*(2d0*lnrtz+5d0/3)
**     print 4,rtz**2*(5d0*lnrtz+59d0/12)
**     print 4,rtz**3*(14d0*lnrtz+449d0/30)
**     print 4,rtz**4*(42d0*lnrtz+1417d0/30)
**     print 4,rtz**5*(132d0*lnrtz+16127d0/105)
**     print 1,bd1ttz
**   &       -1d0/rtz*(b0f(-rtm2,tHmu2,zm2,tm2)+Lnmuzm-1d0)
**     print *
**     print 2,bd2ttz
**     print 4,2d0*lnrtz+5d0/3
**     print 4,rtz*(5d0*lnrtz+59d0/12)
**     print 4,rtz**2*(14d0*lnrtz+449d0/30)
**     print 4,rtz**3*(42d0*lnrtz+1417d0/30)
**     print 4,rtz**4*(132d0*lnrtz+16127d0/105)
**     print 4,rtz**5*(429d0*lnrtz+429697d0/840)
**     print 2,bd2ttz
**   &       -1d0/rtz**2*(b0f_Wdec(-rtm2,tHmu2,zm2,tm2)+Lnmuzm-1d0
**   &                    -rtz*(Lnmutm-Lnmuzm+1d0/2))
**     print *
**     print 3,bdpttz
**     print 4,-lnrtz-17d0/6
**     print 4,-rtz*(7d0*lnrtz+133d0/12)
**     print 4,-rtz**2*(31d0*lnrtz+2657d0/60)
**     print 4,-rtz**3*(126d0*lnrtz+5231d0/30)
**     print 4,-rtz**4*(498d0*lnrtz+142991d0/210)
**     print 4,-rtz**5*(1947d0*lnrtz+2225231d0/840)
**     print 3,bdpttz
**   &       -1d0/rtz*((1d0+2d0*rtz)*rzm2*b0p_Wdec(-rtm2,tm2,zm2)+1d0/2)
*1     format(' bd1ttz=',2f16.12)
*2     format(' bd2ttz=',2f16.12)
*3     format(' bdpttz=',2f16.12)
*4     format(' term  =', f16.12)
**     stop
      ELSE
        bd1ttz=1d0/rtz*(b0f_Wdec(-rtm2,tHmu2,tm2,zm2)+Lnmuzm-1d0)
        bd2ttz=1d0/rtz**2*(b0f_Wdec(-rtm2,tHmu2,tm2,zm2)+Lnmuzm-1d0
     &                    -rtz*(Lnmutm-Lnmuzm+1d0/2))
        bdpttz=1d0/rtz*((1d0+2d0*rtz)*rzm2
     &        *b0p_Wdec(-rtm2,tm2,zm2)+1d0/2)
      ENDIF      
**    print *,'bd1ttz =',bd1ttz
**    print *,'via B0D=',b0d1_Wdec(-rtm2,tHmu2,tm2,zm2)
**    print *,'via B0D=',b0d2_Wdec(-rtm2,tHmu2,tm2,zm2)
**    print *,'bd2ttz =',bd2ttz
**    print *,'via EXT=',1d0/rtz**2*(b0f_Wdec(-rtm2,tHmu2,tm2,zm2)+Lnmuzm-1d0
**   &                   -rtz*(Lnmutm-Lnmuzm+1d0/2))
**    print *,'bdpttz =',bdpttz
**    print *,'via B0D=',b0dp_Wdec(-rtm2,tHmu2,tm2,zm2)
      bd1ttz = b0d1_Wdec(-rtm2,tHmu2,tm2,zm2)
**    bd2ttz = b0d2_Wdec(-rtm2,tHmu2,tm2,zm2)
**    bdpttz = b0dp_Wdec(-rtm2,tHmu2,tm2,zm2)
*-
      IF(rth.LT.1d-3) THEN
        lnrth=DLOG(rth)
        bd1tth=lnrth+1d0/2+rth*(2d0*lnrth+5d0/3)
     &                 +rth**2*(5d0*lnrth+59d0/12)
     &                 +rth**3*(14d0*lnrth+449d0/30)
     &                 +rth**4*(42d0*lnrth+1417d0/30)
     &                 +rth**5*(132d0*lnrth+16127d0/105)
      ELSE
        bd1tth=1d0/rth*(b0f_Wdec(-rtm2,tHmu2,hm2,tm2)+Lnmuhm-1d0)
      ENDIF
**    print *,'bd1tth =',bd1tth
**    print *,'via B0D=',b0d1(-rtm2,tHmu2,tm2,hm2)
      bd1tth = b0d1_Wdec(-rtm2,tHmu2,tm2,hm2)
*-
**    print *,'1,rtw,bw=',rtw,rbw
      IF(rtw.LT.1d-3) THEN
       IF((rbw.GT.1d-20.AND.rbw.LT.1d-3).AND.EWWFFV.EQ.1) THEN
        lnrbw=DLOG(rbw)
        bd1tbw= 1d0/2+rbw+rbw**2+rbw**3+rbw**4
     &         +lnrbw*rbw*(1d0+2d0*rbw+3d0*rbw**2+4d0*rbw**3)
     &  +rtw*(1d0/6+13d0/6*rbw+37d0/6*rbw**2+73d0/6*rbw**3 
     &         +lnrbw*rbw*(1d0+5d0*rbw+14d0*rbw**2))
     &  +rtw**2*(1d0/12+17d0/6*rbw+63d0/4*rbw**2
     &         +lnrbw*rbw*(1d0+9d0*rbw))
     &  +rtw**3*(1d0/20+199d0/60*rbw
     &         +lnrbw*rbw)
     &  +rtw**4*1d0/30
        bdptbw=0d0
*       print 4,bd1tbw
*       test= 1d0/rtw*((1d0-rbw)*(b0f_Wdec(-rtm2,tHmu2,wm2,bm2)+Lnmuwm-1d0)
*    &                     -rbw*(Lnmubm-Lnmuwm))
*       print 5,bd1tbw-test
*5      format(' bd1tbw=',2f15.12)
*       stop
       ELSEIF(rbw.LE.1d-20.OR.EWWFFV.EQ.0) THEN
**      print *,'here'
        bd1tbw=1d0/2+1d0/6*rtw+1d0/12*rtw**2+1d0/20*rtw**3+1d0/30*rtw**4
        bdptbw=-1d0/6-1d0/6*rtw-3d0/20*rtw**2-1d0/6*rtw**3
**      print *,'2,rtw,bw=',rtw,rbw
**      print 4,bd1tbw
**      test= 1d0/rtw*(b0f_Wdec(-rtm2,tHmu2,wm2,bm2)+Lnmuwm-1d0)
**      print 4,bd1tbw-test
**      stop
       ELSE
        PRINT *,'rtw,bw=',rtw,rbw
        PRINT *,'Non-forseen situation with bd1tbw/bdptbw'
        STOP
       ENDIF
      ELSE
       IF(rbw.GT.1d-20) THEN
        bd1tbw=1d0/rtw*((1d0-rbw)*(b0f_Wdec(-rtm2,tHmu2,wm2,bm2)
     &         +Lnmuwm-1d0)
     &                      -rbw*(Lnmubm-Lnmuwm))
        bdptbw=0d0
       ELSEIF(rbw.LE.1d-20.OR.EWWFFV.EQ.0) THEN
        bd1tbw=1d0/rtw*(b0f_Wdec(-rtm2,tHmu2,wm2,bm2)+Lnmuwm-1d0)
        bdptbw=1d0/rtw*(bd1tbw+rwm2*b0p_Wdec(-rtm2,bm2,wm2))
       ENDIF
      ENDIF
**    print *,'tm2,bm2=',tm2,bm2
**    print *,'bd1tbw =',bd1tbw
      bd1tbw = b0dw_Wdec(-rtm2,tHmu2,bm2,wm2)
**    print *,'via B0D=',bd1tbw
 111  CONTINUE
*-
*-  Couplings
*-
      qe =-1d0
      ae =-1d0/2
      i3e=-1d0/2
      ve =-1d0/2-2d0*qe*stw2
      anu=1d0/2
      vnu=1d0/2
*-
      IF    (CHANNL.EQ.0) THEN
        qt = 0d0
        qtc= qt
        qtm= DABS(qt)
        at = 1d0/2
        vt = 1d0/2
        i3t= 1d0/2
        qb =-1d0
        ab =-1d0/2
        vb =-1d0/2-2d0*qb*stw2
        Nc = 1
      ELSEIF(CHANNL.EQ.1.OR.CHANNL.EQ.2.OR.CHANNL.EQ.3) THEN
        qt =-1d0
        qtc= qt
        qtm= DABS(qt)
        at =-1d0/2
        vt =-1d0/2-2d0*qt*stw2
        i3t=-1d0/2
        qb = 0d0
        ab = 1d0/2
        vb = 1d0/2
        Nc = 1
      ELSEIF(CHANNL.EQ.4.OR.CHANNL.EQ.6.OR.CHANNL.EQ.8) THEN
        qt = 2d0/3
        qtc= qt
        qtm= DABS(qt)
        at = 1d0/2
        vt = 1d0/2-2d0*qt*stw2
        i3t= 1d0/2
        qb =-1d0/3
        ab =-1d0/2
        vb =-1d0/2-2d0*qb*stw2
        Nc = 3
      ELSEIF(CHANNL.EQ.5.OR.CHANNL.EQ.7.OR.CHANNL.EQ.9) THEN
        qt =-1d0/3
        qtc= qt
        qtm= DABS(qt)
        at =-1d0/2
        vt =-1d0/2-2d0*qt*stw2
        i3t=-1d0/2
        qb = 2d0/3
        ab = 1d0/2
        vb = 1d0/2-2d0*qb*stw2
        Nc = 3
      ELSEIF(CHANNL.EQ.11.OR.CHANNL.EQ.12.OR.CHANNL.EQ.13) THEN
*-----> qt <--> -qb replacement (see the proof in W_decay_ff_comp.frm)
        qt = 1d0
        qtc= qt
        qtm= DABS(qt)
        at = 1d0/2
        vt = 1d0/2-2d0*qt*stw2
        i3t= 1d0/2
        qb = 0d0
        ab =-1d0/2
        vb =-1d0/2-2d0*qb*stw2
        Nc = 1
      ELSEIF(CHANNL.EQ.14.OR.CHANNL.EQ.15.OR.CHANNL.EQ.16) THEN
        qt = 2d0/3
        qtc= qt
        qtm= DABS(qt)
        at = 1d0/2
        vt = 1d0/2-2d0*qt*stw2
        i3t= 1d0/2
        qb =-1d0/3
        ab =-1d0/2
        vb =-1d0/2-2d0*qb*stw2        
        Nc = 3
      ENDIF
*-
      qe2=qe**2
      qt2=qt**2
      qb2=qb**2
      qeqt=qe*qt
*-
      ae2=ae**2
      ve2=ve**2
*-
      at2=at**2
      vt2=vt**2
*-
      vpae=ve+ae
      vmae=ve-ae
      vpan=1d0
      vman=0d0
      vpat=vt+at
      vmat=vt-at
      vpau=vpat
      vmau=vmat
      vpab=vb+ab
      vmab=vb-ab
      vpad=vpab
      vmad=vmab
      de=vmae/i3e
      dt=vmat/i3t
*-
      IF(FIRST) THEN
        FIRST = .FALSE.
        print *,'Z-mass=',zmA
        print *,'W-mass=',wmA
        print *,'stw**2=',stw2
        print *,'ctw**2=',ctw2
        print *,'H-mass=',hmA
        print *,'alphaI=',alphaI
        print *,'CONHC =',CONHC
        print *,'GammaZ=',GammaZ
      ENDIF
*-
      RETURN 
      END










