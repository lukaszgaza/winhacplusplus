      SUBROUTINE 
     &  BSE_FC_Wdec(Pgg0f,PZgzf,SWWf,Drhof,SWW_fd,SZZ_fd,SHH_fd,Tadfer)
*-------------------------------------------------------------------
*- Fermionic Component of Bosonic Self Energy 
*-------------------
      IMPLICIT NONE!
*-------------------
      INTEGER*4 i,j
      REAL*8 tHmu2
      REAL*8 zm,rwm2,rzm2,rhm2,rmu2,rmd2
      REAL*8 ctw2,stw2,ctw4,ctw6
      REAL*8 rtz,rtw,rhw2,rhz2,rtw2,rhz,rhw,rth,rbw,rbw2,drtwbw
      REAL*8 rmf,rmf2,Qf,cf,I3f,rfw,rfh
      REAL*8 af,af2,vf,vf2
      REAL*8 Lnmuwm,Lnmuzm,Lnmuhm,Lnmutm,Lnmubm
      REAL*8 Pgg0f,PZgzf,SWWf,Drhof,SWW_fd,SZZ_fd,SHH_fd,Tadfer
      REAL*8 b0fwwz,b0fzww,b0fwwh,b0fzzh,b0fww0
      REAL*8 b0pzww,b0pzhz
      REAL*8 b0fhff,b0phff,DREAL,DLOG
      COMPLEX*16 cz2,wm2,zm2,hm2,mf2,mu2,md2
      COMPLEX*16 b0f_Wdec,b0p_Wdec,LOG
      COMPLEX*16 Pgg_f_Wdec,PZg_f_Wdec,SWW_f_Wdec,SZZ_f_Wdec
*-
      COMMON/tHScal/tHmu2
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
      COMMON/BOS_m2/zm,cz2,wm2,zm2,hm2,rwm2,rzm2,rhm2
      COMMON/c_stw2/ctw2,stw2          
      COMMON/c_tw46/ctw4,ctw6             
      COMMON/RATIOS/rtz,rtw,rhw2,rhz2,rtw2,rhz,rhw,rth,rbw,rbw2,drtwbw
      COMMON/SCALES/Lnmuwm,Lnmuzm,Lnmuhm,Lnmutm,Lnmubm
*-
      Pgg0f= DREAL(Pgg_f_Wdec(0d0))
      PZgzf= DREAL(PZg_f_Wdec(rzm2)) 
      SWWf = DREAL(SWW_f_Wdec(rwm2))/rwm2
      Drhof= DREAL(SWW_f_Wdec(rwm2)-SZZ_f_Wdec(rzm2))/rwm2 
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
*---- SE BLOCK FERMIONIC COMPONENT
*-
*---- Tadfer -----------------------------------------------------
      TadFer=0d0
      do 1 i=1,12 
      rmf2=(rmf(i))**2
      TadFer=TadFer-2d0*cf(i)*(DLOG(rmf2/tHmu2)-1d0)*rmf2**2/rwm2/rhm2
 1    continue
*---- SZZ_fd -----------------------------------------------------
      SZZ_fd=0d0
      do 2 i=1,12  
      rmf2=(rmf(i))**2
      af =I3f(i)                       
      vf =I3f(i)-2d0*Qf(i)*stw2   
      af2=af**2                        
      vf2=vf**2  
      mf2=DCMPLX(rmf2,-1d-30)     
      SZZ_fd=SZZ_fd+cf(i)*DREAL(1d0/3*(vf2+af2)*(
     &      -b0f_Wdec(-rzm2,tHmu2,mf2,mf2)+1d0/3
     &      +(2d0*rmf2+rzm2)*b0p_Wdec(-rzm2,mf2,mf2))
     &      -2d0*af2*rmf2*b0p_Wdec(-rzm2,mf2,mf2))
 2    continue
*---- SWW_fd -----------------------------------------------------
      SWW_fd=0d0
      do 3 i=1,6
      rmu2=(rmf(2*i-1))**2
      rmd2=(rmf(2*i)  )**2
      mu2=DCMPLX(rmu2,-1d-30)
      md2=DCMPLX(rmd2,-1d-30)
      if(i.le.3) j=i
      if(i.gt.3) j=i+3
      SWW_fd=SWW_fd+1d0/6d0*cf(j)*DREAL(
     &-(2+(rmu2-rmd2)**2/rwm2**2)*b0f_Wdec(-rwm2,tHmu2,mu2,md2)
     &+(2-(rmu2+rmd2)/rwm2-(rmu2-rmd2)**2/rwm2**2)*rwm2
     &                                         *b0p_Wdec(-rwm2,mu2,md2)
     &+rmu2*(rmd2-rmu2)/rwm2**2*(DREAL(LOG(mu2/tHmu2))-1d0)
     &+rmd2*(rmu2-rmd2)/rwm2**2*(DREAL(LOG(md2/tHmu2))-1d0)
     &+2d0/3)
 3    continue
*---- SHH_fd -----------------------------------------------------
      SHH_fd=0d0
      do 4 i=1,12  
      rmf2=(rmf(i))**2
      mf2=DCMPLX(rmf2,-1d-30) 
      rfw=rmf2/rwm2
      rfh=rmf2/rhm2
      b0fhff = DREAL(b0f_Wdec(-rhm2,tHmu2,mf2,mf2))
      b0phff = DREAL(b0p_Wdec(-rhm2,mf2,mf2))
      SHH_fd=SHH_fd+
     &   cf(i)*rfw/2*(-b0fhff+(1d0-4d0*rfh)*rhm2*b0phff)
 4    continue
*-
      RETURN
      END
