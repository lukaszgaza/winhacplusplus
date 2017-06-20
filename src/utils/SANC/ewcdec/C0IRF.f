      COMPLEX*16 FUNCTION C0IRF_Wdec(M12,M22,P2)
*------------------------------------------
      IMPLICIT NONE!
*-------------------
      REAL*8 P2,D2,DREAL,DABS
      COMPLEX*16 M12,M22,SQR,P2P,P2M,X1,X2,F1
      COMPLEX*16 LOG,SQRT,DCMPLX,XSPENZ_Wdec
      PARAMETER (D2=1.64493406684822618D0)
*-
      IF((DREAL(M12).EQ.DREAL(M22)).AND.DREAL(M12).LT.1d-6) THEN
       C0IRF_Wdec=1D0/P2
     &           *(1D0/2D0*(LOG((P2-DCMPLX(0D0,1D-20))/M12))**2-D2)
       RETURN
      ENDIF
*-
      SQR=SQRT(P2**2+2D0*P2*(M12+M22)+(M12-M22)**2)
        IF(DABS(DREAL(P2P)).GE.DABS(DREAL(P2M))) THEN
        P2P=P2-M12+M22+SQR
        P2M=P2-M12+M22-SQR
        IF(DREAL(P2P).GE.DREAL(P2M)) THEN
         X1=P2P/2D0/P2
         X2=-2D0*M12/P2P
        ELSE
         X2=P2M/2D0/P2
         X1=-2D0*M12/P2M
        ENDIF
      ELSE
        P2P=P2-M22+M12+SQR
        P2M=P2-M22+M12-SQR
        IF(DABS(DREAL(P2P)).GE.DABS(DREAL(P2M))) THEN
         X1=P2P/2D0/P2
         X2=-2D0*M22/P2P
        ELSE
         X2=P2M/2D0/P2
         X1=-2D0*M22/P2M
        ENDIF
      ENDIF
*-
      F1=1D0/SQR*(LOG(1D0-1D0/X2)-LOG(1D0-1D0/X1))
*-
      C0IRF_Wdec=1D0/2*F1*LOG((P2-DCMPLX(0D0,1D-20))/M22)
     &    +1D0/2*(1D0/SQR*(
     & +.5D0*LOG((1D0-X2)*(X1-X2))**2
     & -.5D0*LOG(   (-X2)*(X1-X2))**2
     & -.5D0*LOG((X1-1D0)*(X1-X2))**2
     & +.5D0*LOG( X1     *(X1-X2))**2
     & -XSPENZ_Wdec((1D0-X2)/(X1-X2))+XSPENZ_Wdec((-X2)/(X1-X2))
     & +XSPENZ_Wdec((1D0-X1)/(X2-X1))-XSPENZ_Wdec((-X1)/(X2-X1))
     &                    ))
*-
      RETURN    
      END
