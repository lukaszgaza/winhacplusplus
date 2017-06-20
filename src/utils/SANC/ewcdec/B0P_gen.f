      COMPLEX*16 FUNCTION B0P_Wdec(Q2,M12,M22)
*     -----------------------------------
* general B0P; CalcPHEP group
*
      IMPLICIT NONE
      INTEGER*4 N
      REAL*8 Q2,ABQ2,DRM12,DRM22,RU,RD,B0PEXP_Wdec,DREAL,DABS
      COMPLEX*16 M12,M22,SD,B,XPL,XMI,FB0P_Wdec,SUM,DCMPLX,LOG,SQRT
*
      ABQ2=DABS(Q2)
      DRM12=DREAL(M12)
      DRM22=DREAL(M22)
      SUM =DCMPLX(0D0,0D0)
      IF(ABQ2.NE.0D0) THEN
* B0P(Q2;M12,M22)
       SD=SQRT(Q2**2+2D0*Q2*(M12+M22)+(M12-M22)**2)
       B=Q2-M22+M12
       XPL=(-B+SD)/(-2D0*Q2)
       IF(DABS(DREAL(XPL)).LT.1D-7) XPL=-2D0*M22/(B+SD)
       XMI=(-B-SD)/(-2D0*Q2)
       IF(DABS(DREAL(XMI)).LT.1D-7) XMI=-2D0*M22/(B-SD)
       IF(DRM12.GT.1D-14.AND.DRM22.GT.1D-14) THEN
        IF(ABQ2.LE.1D0.OR.DRM12.LE.1D0.OR.DRM22.LE.1D0) THEN 
* TAYLOR EXPANSIONS
         IF((ABQ2/DRM12.LT.1D-3).AND.(DRM22/DRM12.LT.1D-3)) THEN
          RU=ABQ2/DRM12
          RD=DRM22/DRM12
          B0P_Wdec=1D0/DRM12*B0PEXP_Wdec(RU,RD)
         ELSEIF((ABQ2/DRM22.LT.1D-3).AND.(DRM12/DRM22.LT.1D-3)) THEN
          RU=ABQ2/DRM22
          RD=DRM12/DRM22
          B0P_Wdec=1D0/DRM22*B0PEXP_Wdec(RU,RD)
         ELSE
          DO N=1,2
           SUM=SUM+(-1)**N*(FB0P_Wdec(N,XPL)-FB0P_Wdec(N,XMI))
          ENDDO
          B0P_Wdec=-1D0/SD*SUM
         ENDIF
        ELSE
         DO N=1,2
          SUM=SUM+(-1)**N*(FB0P_Wdec(N,XPL)-FB0P_Wdec(N,XMI))
         ENDDO
         B0P_Wdec=-1D0/SD*SUM
        ENDIF
       ELSEIF(DRM12.LE.1D-14.AND.DRM22.GT.1D-14) THEN
        IF(ABQ2.LE.1D0) THEN
         RU=ABQ2/DRM22
         B0P_Wdec=-1D0/DRM22*(1D0/2+1D0/3*RU+1D0/4*RU**2+1D0/5*RU**3
     &                  +1D0/6*RU**4+1D0/7*RU**5+1D0/8*RU**6)
        ELSE
         B0P_Wdec=-1D0/Q2**2*(Q2-M22*LOG(1D0+Q2/M22))
        ENDIF 
      ELSEIF(DRM12.GT.1D-14.AND.DRM22.LE.1D-14) THEN
        IF(ABQ2.LE.1D0) THEN
         RU=ABQ2/DRM12
         B0P_Wdec=-1D0/DRM12*(1D0/2+1D0/3*RU+1D0/4*RU**2+1D0/5*RU**3
     &                  +1D0/6*RU**4+1D0/7*RU**5+1D0/8*RU**6)
        ELSE
         B0P_Wdec=-1D0/Q2**2*(Q2-M12*LOG(1D0+Q2/M12))
        ENDIF      
       ELSEIF(DRM12.LE.1D-14.AND.DRM22.LE.1D-14) THEN
        B0P_Wdec=-1D0/Q2
       ELSE
        PRINT *,'NOT FORESEEN SET OF MASSES IN B0P: M12,M22=',M12,M22
        STOP
       ENDIF        
      ELSE
* B0P(0;M12,M22)
       IF((DRM12.NE.DRM22).AND.DRM22.NE.0D0.AND.DRM12.NE.0D0) THEN
        B0P_Wdec=M12*M22/(M12-M22)**3*LOG(M12/M22)
     &              -(M12+M22)/2D0/(M12-M22)**2
       ELSEIF(DRM12.EQ.0D0.AND.DRM22.NE.0D0) THEN
        B0P_Wdec=-1D0/2D0/M22
       ELSEIF(DRM12.NE.0D0.AND.DRM22.EQ.0D0) THEN
        B0P_Wdec=-1D0/2D0/M12
       ELSEIF(DRM12.EQ.DRM22) THEN
        B0P_Wdec=-1D0/6D0/M12
       ELSE
        PRINT*,'B0P(0;...) NOT FORESEEN SET OF MASSES: M12,M22=',M12,M22
        STOP        
       ENDIF
      ENDIF
*
      RETURN    
      END

      COMPLEX*16 FUNCTION FB0P_Wdec(N,Y)
*
      IMPLICIT NONE
      INTEGER*4 N,L
      COMPLEX*16 Y,SUM,DCMPLX
*
      SUM=DCMPLX(0D0,0D0)
      DO L=1,N
        SUM=SUM+Y**(N-L)/L
      ENDDO
      FB0P_Wdec=-Y**N*LOG(1D0-1D0/Y)-SUM
*
      RETURN
      END

      REAL*8 FUNCTION B0PEXP_Wdec(RU,RD)
*
      IMPLICIT NONE
      REAL*8 RU,RD,LNRD
*
      LNRD=DLOG(RD)
      B0PEXP_Wdec=
     &-1D0/2*(1D0        +3D0*RD+  5D0*RD**2+    7D0*RD**3+   9D0*RD**4
     &     +11D0*RD**5)
     &     -LNRD*RD*(1D0+ 3D0*RD+  6D0*RD**2+   10D0*RD**3+  15D0*RD**4)
     &+RU   *(-1D0/3-  14D0/3*RD- 17D0*RD**2-124D0/3*RD**3-245D0/3*RD**4 
     &     -LNRD*RD*(2D0+12D0*RD+ 40D0*RD**2+  100D0*RD**3) )
     &+RU**2*(-1D0/4  -35D0/4*RD- 56D0*RD**2-  210D0*RD**3 
     &     -LNRD*RD*(3D0+30D0*RD+150D0*RD**2) )
     &+RU**3*(-1D0/5-202D0/15*RD-134D0*RD**2 
     &     -LNRD*RD*(4D0+60D0*RD) )
     &+RU**4*(-1D0/6-56D0/3*RD-5D0*LNRD*RD)
     &+RU**5*(-1D0/7)  
*
      RETURN
      END
