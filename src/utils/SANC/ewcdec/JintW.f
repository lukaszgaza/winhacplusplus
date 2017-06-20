      COMPLEX*16 FUNCTION JINTW_Wdec(M12,M22,P2)
*------------------------------------------
      IMPLICIT NONE!
*-------------------
      REAL*8 P2,DABS
      COMPLEX*16 M12,M22,SQR,P2P,P2M,X1,X2,SQRT,LOG
*-
      SQR=SQRT(P2**2+2D0*P2*(M12+M22)+(M12-M22)**2)
***   X1=(P2-M12+M22+SQR)/2D0/P2
***   X2=(P2-M12+M22-SQR)/2D0/P2
      IF(DREAL(M12).LE.DREAL(M22)) THEN
        P2P=P2-M12+M22+SQR
        P2M=P2-M12+M22-SQR
        IF(DABS(DREAL(P2P)).GE.DABS(DREAL(P2M))) THEN
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
      JINTW_Wdec=-(LOG(1D0-1D0/X2)-LOG(1D0-1D0/X1))
*-
      RETURN
      END
