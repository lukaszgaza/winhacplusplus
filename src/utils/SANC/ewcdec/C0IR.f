      COMPLEX*16 FUNCTION C0IR_Wdec(tHla2,M12,M22,P2)
*     ------------------------------------------
* All cases of C0IR: 
*
      IMPLICIT NONE
      INTEGER*4 Lim
      REAL*8 tHla2,P2,DRM12,DRM22,D2,DREAL
      COMPLEX*16 M12,M22,SQR,P2P,P2M,X1,X2,OMX1,OMX2,F1
      COMPLEX*16 LOG,SQRT,DCMPLX,XSPENZ_Wdec
      PARAMETER (D2=1.64493406684822618D0)
*
      DRM12=DREAL(M12)
      DRM22=DREAL(M22) 
*
* Limiting cases
*
      IF(DRM12.LE.1D-2.AND.DRM22.LE.1D-2) THEN
* M1.EQ.M2
       IF(DRM12.EQ.DRM22) THEN
        C0IR_Wdec=1D0/P2*(LOG(M12/tHla2)*LOG((P2-DCMPLX(0D0,1D-20))/M12)
     &      +1D0/2D0*(LOG((P2-DCMPLX(0D0,1D-20))/M12))**2-D2)
        RETURN
* M1.NE.M2
       ELSE
        C0IR_Wdec=1D0/2/P2*(
     &      +(LOG(P2/M12)+LOG(P2/M22))*LOG((P2-DCMPLX(0D0,1D-20))/tHla2)
     &      -.5D0*LOG(M12/P2)**2-.5D0*LOG(M22/P2)**2-2D0*D2)
        RETURN
       ENDIF
      ENDIF
*
* Special case for W->\nu\lept decay
*
      IF(DABS(P2).LT.1D-10) THEN
        C0IR_Wdec=1D0/4/DREAL(M12-M22)*DLOG(DRM12/DRM22)
     &       *DLOG(DRM12*DRM22/tHla2**2)
        RETURN
      ENDIF
*
* General, exact case of C0IR
*
      SQR=SQRT(P2**2+2D0*P2*(M12+M22)+(M12-M22)**2)
      P2P=P2-M12+M22+SQR
      P2M=P2-M12+M22-SQR
      IF(DABS(DREAL(P2P)).GE.DABS(DREAL(P2M))) THEN
        X1=P2P/2D0/P2
        X2=-2D0*M12/P2P
      ELSE
        X2=P2M/2D0/P2
        X1=-2D0*M12/P2M
      ENDIF
      P2P=P2+M12-M22+SQR  
      P2M=P2+M12-M22-SQR
      IF(DABS(DREAL(P2P)).GE.DABS(DREAL(P2M))) THEN
        OMX2=P2P/2D0/P2
        OMX1=-2D0*M22/P2P
      ELSE
        OMX1=P2M/2D0/P2
        OMX2=-2D0*M22/P2M
      ENDIF
*
      F1=1D0/SQR*(LOG(OMX2/(-X2))-LOG(OMX1/(-X1)))
*
      C0IR_Wdec=1D0/2*F1*LOG((P2-DCMPLX(0D0,1D-20))/tHla2)  
     &    +1D0/2/SQR*(
     & +.5D0*LOG( OMX2*(X1-X2))**2
     & -.5D0*LOG(  -X2*(X1-X2))**2
     & -.5D0*LOG(-OMX1*(X1-X2))**2
     & +.5D0*LOG(   X1*(X1-X2))**2
     & -XSPENZ_Wdec(OMX2/(X1-X2))+XSPENZ_Wdec((-X2)/(X1-X2))
     & +XSPENZ_Wdec(OMX1/(X2-X1))-XSPENZ_Wdec((-X1)/(X2-X1))
     &               )
*
      RETURN    
      END
