      COMPLEX*16 FUNCTION BFF_Wdec(Q2,MU2,M12,M22)
*     ---------------------------------------
* general BFF; 
*
      IMPLICIT NONE
      REAL*8 MU2,Q2,ABQ2,DRM12,DRM22,RU,RD,BFFEXP_Wdec,DREAL,DABS,DLOG
      COMPLEX*16 M12,M22,B,SQR,SQS,LOG,SQRT,B0F_Wdec
*
      ABQ2=DABS(Q2)
      SQR=SQRT(Q2**2+2D0*Q2*(M12+M22)+(M12-M22)**2)
      SQS=Q2+M12+M22
      DRM12=DREAL(M12)
      DRM22=DREAL(M22)
*
      IF(ABQ2.NE.0D0) THEN
*
* BFF(Q2;MU2,M12,0)
       IF(DRM22.LT.1D-14) THEN
        BFF_Wdec=1D0/3D0*LOG(M12/MU2)
     &          -5D0/9D0+2D0/3D0*(M12/Q2+M12**2/Q2**2)
     &   +(1D0-M12/Q2-2D0*M12**2/Q2**2)/3D0*(1D0+M12/Q2)*LOG(1D0+Q2/M12)
        RETURN
       ENDIF
*
* BFF(Q2;MU2,0,M22)
       IF(DRM12.LT.1D-14) THEN
         BFF_Wdec=1D0/3D0*LOG(M22/MU2)
     &           -5D0/9D0+2D0/3D0*(M22/Q2+M22**2/Q2**2)
     &   +(1D0-M22/Q2-2D0*M22**2/Q2**2)/3D0*(1D0+M22/Q2)*LOG(1D0+Q2/M22)
        RETURN
       ENDIF
*
* TAYLOR EXPANSIONS: TO STABILIZE BFF(-SMALL,MU2,SMALL,M22)
       IF((ABQ2/DRM12.LT.1D-3).AND.(DRM22/DRM12.LT.1D-3)) THEN
        RU=(-Q2)/DRM12
        RD=DRM22/DRM12
        BFF_Wdec=BFFEXP_Wdec(RU,RD)+1D0/3*DLOG(DRM12/MU2)
        RETURN
       ELSEIF((ABQ2/DRM22.LT.1D-3).AND.(DRM12/DRM22.LT.1D-3)) THEN
        RU=(-Q2)/DRM22
        RD=DRM12/DRM22
        BFF_Wdec=BFFEXP_Wdec(RU,RD)+1D0/3*DLOG(DRM22/MU2)
        RETURN
       ENDIF
*
* BFF(Q2;MU2,M12,M22)
       IF(DRM12.GT.1D-14.AND.DRM22.GT.1D-14) THEN
          BFF_Wdec=1d0/3/Q2*
     &   (Q2/3-2D0*(M12-M22)**2/Q2
     &   +(1D0+2D0*(M12-M22)/Q2)*M12*LOG(M12/MU2)
     &   +(1D0+2D0*(M22-M12)/Q2)*M22*LOG(M22/MU2)
     &   -(Q2-M12-M22-2D0*(M12-M22)**2/Q2)*B0F_Wdec(Q2,MU2,M12,M22))  
       ELSEIF(DRM12.LE.1D-14.AND.DRM22.GT.1D-14) THEN
         BFF_Wdec=1D0/3D0*LOG(M22/MU2)
     &           -5D0/9D0+2D0/3D0*(M22/Q2+M22**2/Q2**2)
     &   +(1D0-M22/Q2-2D0*M22**2/Q2**2)/3D0*(1D0+M22/Q2)*LOG(1D0+Q2/M22)
       ELSEIF(DRM12.GT.1D-14.AND.DRM22.LE.1D-14) THEN
         BFF_Wdec=1D0/3D0*LOG(M12/MU2)
     &           -5D0/9D0+2D0/3D0*(M12/Q2+M12**2/Q2**2)
     &   +(1D0-M12/Q2-2D0*M12**2/Q2**2)/3D0*(1D0+M12/Q2)*LOG(1D0+Q2/M12)
       ELSEIF(DRM12.LE.1D-14.AND.DRM22.LE.1D-14) THEN
         BFF_Wdec=1D0/3D0*LOG(DCMPLX(Q2,-1D-20)/MU2)-5D0/9D0
       ELSE
        PRINT *,'NOT FORESEEN SET OF MASSES IN BFF: M12,M22=',M12,M22
        STOP
       ENDIF
*
      ELSEIF(ABQ2.EQ.0D0) THEN
*
* BFF(0;MU2,M12,M22)
       IF((DRM12.NE.DRM22).AND.DRM22.NE.0D0.AND.DRM12.NE.0D0) THEN
        B=M12/M22-1D0
        BFF_Wdec=1D0/3D0*LOG(M22/MU2)
     &     +(1D0/3D0-1D0/B**2-2D0/3D0/B**3)*LOG(M12/M22)
     &     -5D0/18D0+2D0/3D0/B+2D0/3D0/B**2
       ELSEIF(DRM12.EQ.0D0.AND.DRM22.NE.0D0) THEN
        BFF_Wdec=1D0/3D0*LOG(M22/MU2)-5D0/18D0
       ELSEIF(DRM12.NE.0D0.AND.DRM22.EQ.0D0) THEN
        BFF_Wdec=1D0/3D0*LOG(M12/MU2)-5D0/18D0
       ELSEIF(DRM12.EQ.DRM22) THEN
        BFF_Wdec=1D0/3D0*LOG(M12/MU2)
       ELSE
        PRINT*,'BFF(0;...) NOT FORESEEN SET OF MASSES: M12,M22=',M12,M22
        STOP        
       ENDIF
      ELSE
       PRINT*,'B0F(Q2;...) NOT FORESEEN SET OF MASSES: M12,M22=',M12,M22
      ENDIF           
*
      RETURN    
      END  

      REAL*8 FUNCTION BFFEXP_Wdec(RU,RD)
*
      IMPLICIT NONE
      REAL*8 RU,RD,LNRD
*
      LNRD=DLOG(RD)
*
      BFFEXP_Wdec=
     &  - 5D0/18 + 2D0/3*RD + 4D0/3*RD**2 + 2D0*RD**3 + 8D0/3*RD**4 
     &  + 1760D0/3*RD**5 
     &  + LNRD*RD**2*(1D0+8D0/3*RD+5D0*RD**2+792D0*RD**3)
     &  + RU*(-1D0/6+1D0/2*RD+25D0/6*RD**2+77D0/6*RD**3+17074D0/9*RD**4
     &  + LNRD*RD**2*(2D0+10D0*RD+4990D0/3*RD**2))
     &  + RU**2*(-1D0/20+2D0/5*RD+157D0/20*RD**2+98506D0/45*RD**3
     &  + LNRD*RD**2*(3D0+3992D0/3*RD))
     &  + RU**3*(-1D0/45+1D0/3*RD+13769D0/15*RD**2+396D0*LNRD*RD**2)
     &  + RU**4*(-1D0/84+4336D0/35*RD+112D0/3*LNRD*RD)
     &  + RU**5*(53D0/15+2D0/3*LNRD)
*
      RETURN
      END


    
