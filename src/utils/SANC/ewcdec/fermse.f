      SUBROUTINE FERMSE_Wdec(I,s,Pggf,PZgf,DZf,DWf)   
*---------------------------------------------
      IMPLICIT NONE!
*-------------------
      INTEGER*4 I
      REAL*8 s,zm,rwm2,rzm2,rhm2
      REAL*8 ctw2,stw2,DREAL
      COMPLEX*16 cz2,wm2,zm2,hm2
      COMPLEX*16 Pggf,PZgf,DZf,DWf
      COMPLEX*16 Pgg_f_Wdec,PZg_f_Wdec,DW_f_Wdec,SZZ_f_Wdec,DCMPLX
*-
      COMMON/BOS_m2/zm,cz2,wm2,zm2,hm2,rwm2,rzm2,rhm2
      COMMON/c_stw2/ctw2,stw2
*-
      Pggf = Pgg_f_Wdec(s)
      PZgf = PZg_f_Wdec(s)
      DZf  = DCMPLX(0d0,0d0)
      IF(I.EQ.1) THEN
       IF(rzm2-s.EQ.0d0) THEN
        DZf = (SZZ_f_Wdec(s*(1d0+1d-10))-SZZ_f_Wdec(rzm2))/ctw2/
     &        (rzm2-s*(1d0+1d-10))
       ELSE
        DZf = (SZZ_f_Wdec(s)-SZZ_f_Wdec(rzm2))/ctw2/(rzm2-s)
       ENDIF
      ENDIF
      DWf=DW_f_Wdec(s)
*-    
      RETURN
      END

      COMPLEX*16 FUNCTION DW_f_Wdec(s)
*--------------------------------
      IMPLICIT NONE!
*-------------------
      INTEGER*4 i
      REAL*8 tHmu2
      REAL*8 s,zm,rwm2,rzm2,rhm2
      REAL*8 rmf,Qf,cf,I3f
      REAL*8 DLOG,DREAL
      REAL*8 bffwdu,b1fwud,b1fwdu,b0fwdu
      COMPLEX*16 bff_Wdec,bffsdu,b1fsud,b1fsdu,b0f_Wdec,b0fsdu,b0p_Wdec
      COMPLEX*16 cz2,wm2,zm2,hm2,mu2,md2,DCMPLX
*-
      COMMON/tHScal/tHmu2
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
      COMMON/BOS_m2/zm,cz2,wm2,zm2,hm2,rwm2,rzm2,rhm2
*-
      DW_f_Wdec=DCMPLX(0d0,0d0)
      DO 1 i=1,6
      mu2=DCMPLX((rmf(2*i-1))**2,-1d-20)
      md2=DCMPLX((rmf(2*i)  )**2,-1d-20)
      bffwdu=DREAL(bff_Wdec(-rwm2,tHmu2,md2,mu2))
      b0fwdu=DREAL(b0f_Wdec(-rwm2,tHmu2,md2,mu2))
      bffsdu=bff_Wdec(-s,tHmu2,md2,mu2)
      b0fsdu=b0f_Wdec(-s,tHmu2,md2,mu2)
      b1fwdu=DREAL(-1d0/2/rwm2*         
     &        (md2*(DLOG(DREAL(md2)/tHmu2)-1d0)                
     &        -mu2*(DLOG(DREAL(mu2)/tHmu2)-1d0)        
     &        +(md2-mu2+rwm2)*b0fwdu))
      b1fwud=DREAL(-1d0/2/rwm2*
     &        (mu2*(DLOG(DREAL(mu2)/tHmu2)-1d0)                
     &        -md2*(DLOG(DREAL(md2)/tHmu2)-1d0)       
     &        +(mu2-md2+rwm2)*b0fwdu)) 
      IF(s.NE.0d0) THEN
       b1fsdu=-1d0/2/s*            
     &        (md2*(DLOG(DREAL(md2)/tHmu2)-1d0)                
     &        -mu2*(DLOG(DREAL(mu2)/tHmu2)-1d0)        
     &        +(md2-mu2+s)*b0fsdu)
       b1fsud=-1d0/2/s*
     &        (mu2*(DLOG(DREAL(mu2)/tHmu2)-1d0)                
     &        -md2*(DLOG(DREAL(md2)/tHmu2)-1d0)       
     &        +(mu2-md2+s)*b0fsdu) 
      ELSE
       b1fsdu=-1d0/2*b0f_Wdec(0d0,tHmu2,md2,mu2)
     &        +1d0/2*(md2-mu2)*b0p_Wdec(0d0,md2,mu2)
       b1fsud=-1d0/2*b0f_Wdec(0d0,tHmu2,mu2,md2)
     &        +1d0/2*(mu2-md2)*b0p_Wdec(0d0,mu2,md2)
      ENDIF        
      DW_f_Wdec =DW_f_Wdec +( - cf(2*i)*(s*bffsdu-rwm2*bffwdu)
     &         +cf(2*i-1)*mu2*(  b1fsdu-     b1fwdu)
     &         +cf(2*i)  *md2*(  b1fsud-     b1fwud))/(-s+rwm2)
 1    CONTINUE
      END

      COMPLEX*16 FUNCTION SWW_f_Wdec(s)
*---------------------------------
      IMPLICIT NONE!
*-------------------
      INTEGER*4 i
      REAL*8 s,tHmu2
      REAL*8 rmf,Qf,cf,I3f
      REAL*8 DLOG
      COMPLEX*16 bff_Wdec,bffsdu,b1fsud,b1fsdu,b0f_Wdec,b0fsdu,b0p_Wdec
      COMPLEX*16 mu2,md2,DCMPLX
*-
      COMMON/tHScal/tHmu2
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
*-
      SWW_f_Wdec=DCMPLX(0d0,0d0)
      DO 1 i=1,6
      mu2=DCMPLX((rmf(2*i-1))**2,-1d-20)
      md2=DCMPLX((rmf(2*i)  )**2,-1d-20)
      bffsdu=bff_Wdec(-s,tHmu2,md2,mu2)
      b0fsdu=b0f_Wdec(-s,tHmu2,md2,mu2)
      IF(s.NE.0d0) THEN
       b1fsdu=-1d0/2/s*            
     &        (md2*(DLOG(DREAL(md2)/tHmu2)-1d0)                
     &        -mu2*(DLOG(DREAL(mu2)/tHmu2)-1d0)        
     &        +(md2-mu2+s)*b0fsdu)
       b1fsud=-1d0/2/s*
     &        (mu2*(DLOG(DREAL(mu2)/tHmu2)-1d0)                
     &        -md2*(DLOG(DREAL(md2)/tHmu2)-1d0)       
     &        +(mu2-md2+s)*b0fsdu) 
      ELSE
       b1fsdu=-1d0/2*b0f_Wdec(0d0,tHmu2,md2,mu2)
     &        +1d0/2*(md2-mu2)*b0p_Wdec(0d0,md2,mu2)
       b1fsud=-1d0/2*b0f_Wdec(0d0,tHmu2,mu2,md2)
     &        +1d0/2*(mu2-md2)*b0p_Wdec(0d0,mu2,md2)
      ENDIF        
      SWW_f_Wdec =SWW_f_Wdec - s*cf(2*i)*bffsdu
     &       +cf(2*i-1)*mu2*b1fsdu+cf(2*i)*md2*b1fsud
 1    CONTINUE
      END

      COMPLEX*16 FUNCTION SZZ_f_Wdec(s)
*---------------------------------
      IMPLICIT NONE!
*-------------------
      INTEGER*4 i
      REAL*8 s,tHmu2
      REAL*8 rmf,rmf2,Qf,cf,I3f
      REAL*8 af,af2,vf,vf2
      REAL*8 ctw2,stw2
      COMPLEX*16 bff_Wdec,b0f_Wdec,bffsff,b0fsff,mf2,DCMPLX
*-
      COMMON/tHScal/tHmu2
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
      COMMON/c_stw2/ctw2,stw2
*-
      SZZ_f_Wdec=DCMPLX(0d0,0d0)
      DO 1 i=1,12
      rmf2=(rmf(i))**2
      mf2 =DCMPLX(rmf2,-1d-20)     
      bffsff = bff_Wdec(-s,tHmu2,mf2,mf2)   
      b0fsff = b0f_Wdec(-s,tHmu2,mf2,mf2)   
      af =I3f(i)                       
      vf =I3f(i)-2d0*Qf(i)*stw2   
      af2=af**2                        
      vf2=vf**2
      SZZ_f_Wdec=SZZ_f_Wdec+cf(i)*(-(vf2+af2)*s*bffsff
     &          -2d0*af2*rmf2*b0fsff)
 1    CONTINUE
      END

      COMPLEX*16 FUNCTION SHH_f_Wdec(s)
*---------------------------------
      IMPLICIT NONE!
*-------------------
      INTEGER*4 i
      REAL*8 tHmu2
      REAL*8 s,zm,rwm2,rzm2,rhm2
      REAL*8 rmf,rmf2,Qf,cf,I3f
      COMPLEX*16 cz2,wm2,zm2,hm2
      COMPLEX*16 b0f_Wdec,b0fsff,mf2,DCMPLX
*-
      COMMON/tHScal/tHmu2
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
*-
      COMMON/BOS_m2/zm,cz2,wm2,zm2,hm2,rwm2,rzm2,rhm2
*-
      SHH_f_Wdec=DCMPLX(0d0,0d0)
      DO 1 i=1,12
      rmf2=(rmf(i))**2
      mf2 =DCMPLX(rmf2,-1d-20)     
      b0fsff = b0f_Wdec(-s,tHmu2,mf2,mf2)   
      SHH_f_Wdec=SHH_f_Wdec+cf(i)*rmf2/rwm2*(!rmf2*(LOG(mf2/tHmu2)-1)
     &           +(s-4d0*rmf2)/2d0*b0fsff)
 1    CONTINUE
      END

      COMPLEX*16 FUNCTION Pgg_f_Wdec(s)
*---------------------------------
      IMPLICIT NONE!
*-------------------
      INTEGER*4 i
      REAL*8 s,tHmu2
      REAL*8 rmf,rmf2,Qf,cf,I3f
      COMPLEX*16 bff_Wdec,bffsff,mf2,DCMPLX
*-
      COMMON/tHScal/tHmu2
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
*-
      Pgg_f_Wdec =DCMPLX(0d0,0d0)
      DO 1 i=1,12  
      rmf2=(rmf(i))**2
      mf2=DCMPLX(rmf2,-1d-20)     
      bffsff=bff_Wdec(-s,tHmu2,mf2,mf2)
      Pgg_f_Wdec =Pgg_f_Wdec+4d0*cf(i)*Qf(i)**2*bffsff
 1    CONTINUE
      END

      COMPLEX*16 FUNCTION PZg_f_Wdec(s)
*---------------------------------
      IMPLICIT NONE!
*-------------------
      INTEGER*4 i
      REAL*8 s,tHmu2
      REAL*8 rmf,Qf,Qfm,cf,I3f
      REAL*8 ctw2,stw2 
      COMPLEX*16 bff_Wdec,bffsff,mf2,DCMPLX
*-
      COMMON/tHScal/tHmu2
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
      COMMON/c_stw2/ctw2,stw2
*-
      PZg_f_Wdec =DCMPLX(0d0,0d0)
      DO 1 i=1,12
      mf2=DCMPLX((rmf(i))**2,-1d-20)     
      bffsff=bff_Wdec(-s,tHmu2,mf2,mf2)   
      Qfm   =abs(Qf(i))
      PZg_f_Wdec =PZg_f_Wdec+cf(i)*Qfm*(1d0-4d0*stw2*Qfm)*bffsff
 1    CONTINUE
      END
