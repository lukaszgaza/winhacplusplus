      COMPLEX*16 FUNCTION C0_Wdec(p12,p22,Q2,m12,m22,m32)
*     ----------------------------------------------
* general C0; CalpHEP group, version 17/01/01
* Attention!  It doesn't work under physical treshold! To be improved!!! 
*
      IMPLICIT NONE
      INTEGER I
      REAL*8 eps,p12,p22,Q2,p1dp2,DET3,SQD3,alpha,N,d,e,r0,pi,zet2
      REAL*8 rm12,rm22,rm32,rmb2,YR,RLOG,TMI,YPL,YMI,YDE,SIG
      REAL*8 DABS,DSQRT,DLOG,DDILOG_Wdec,DREAL,DIMAG
      COMPLEX*16 m12,m22,m32,mb2,ratio,LAMBDA,ALS,BETA
      COMPLEX*16 SQS,SQR,X0,X1,X2,X3,Y0,Y1,Y2,Y3,Y4,Y5,Y6,XSPENZ_Wdec
      COMPLEX*16 R21,R31,R13,R23,BETA21,BETA31,BETA13,BETA23,BETAZ
      COMPLEX*16 IMI,RKERN_Wdec,SKERN,SKERN1,SKERN2
      COMPLEX*16 AP,BP,CP,DP,SQDP,YP,YM,EP,DE,TERM,VELTET_Wdec
      COMPLEX*16 SQRT,LOG,DCMPLX
*
      DIMENSION AP(0:2),BP(0:2),CP(0:2),DP(0:2),SQDP(0:2)
      DIMENSION EP(0:2),DE(0:2),YR(0:2),YP(0:2),YM(0:2)
      DIMENSION SKERN1(0:2),SKERN2(0:2),SKERN(0:2)
*
      PARAMETER (eps=1D-20)
      PARAMETER (pi=3.14159265358979D0)
      PARAMETER (zet2=1.64493406684822618D0)
*
      rm12=DREAL(m12)
      rm22=DREAL(m22)
      rm32=DREAL(m32)
*
* Special cases, moved from C0_sp
*
*     c0(0,-rem2,Qs;wm2,em2,0)
*     c0(-rem2,0,Qs;0,wm2,em2)
*     c0(-rum2,0,Qs;0,um2,wm2) ???
*
*     print *,'args=',p12,p22,Q2,m12,m22,m32
      IF((p12.EQ.0D0.AND.-p22.EQ.rm22.AND.rm32.EQ.0D0).OR.
     &   (-p12.EQ.rm22.AND.p22.EQ.0D0.AND.rm12.EQ.0D0)) THEN
       IF(p12.EQ.0D0.AND.rm32.EQ.0D0) THEN
        mb2=m12
       ELSEIF(p22.EQ.0D0.AND.rm12.EQ.0D0) THEN
        mb2=m32
       ENDIF
       IF(rm22.LT.1D-2) THEN
        ratio=mb2/(Q2+mb2)
        C0_Wdec=1D0/Q2*(-1D0/2*LOG(ratio)**2
     &             -LOG(Q2/m22)*LOG(ratio)-XSPENZ_Wdec(ratio)+zet2)
        RETURN
       ELSE
        C0_Wdec=1D0/(Q2+m22)*(XSPENZ_Wdec(-Q2/m22*(mb2-m22)/(Q2+mb2))
     &    -XSPENZ_Wdec(-Q2/m22)-XSPENZ_Wdec((mb2-m22)/(Q2+mb2))+zet2)
        RETURN
       ENDIF
*
*     c0(-rem2,-rem2,Qs;cz2,em2,cz2)
*     c0(-rem2,-rem2,Qs;zm2,em2,cz2)
*     c0(-rem2,-rem2,Qs;cz2,em2,zm2)
*
      ELSEIF(p12.EQ.p22.AND.-p12.EQ.rm22.AND.(rm12.EQ.0D0.OR.rm32.EQ.0D0
     &                                       )) THEN
       IF(rm12.NE.0D0) mb2=m12
       IF(rm32.NE.0D0) mb2=m32
       IF(rm12.EQ.0D0.AND.rm32.EQ.0D0) mb2=m32
       rmb2=DREAL(mb2)
       IF(rmb2.LT.1D-2.AND.rm22.LT.1D-2) THEN
* GG light mass case ( LOG(Q2/M2)=LOG(S/M2)+SIG*DCMPLX(0D0,1D0)*PI )
        SIG=DABS(Q2)/Q2
        ALS=LOG(-Q2/m22)
        C0_Wdec=1D0/Q2*(1D0/2*ALS**2+SIG*DCMPLX(0D0,1D0)*PI*ALS+zet2)
        RETURN
       ELSEIF(rmb2.LT.1D-2.AND.rm22.GE.1D-2) THEN
* GG massive case
        BETA=SQRT(1D0+4D0*m22/Q2)
        Y1=-Q2/(2D0*m22)*(1D0+BETA)
        Y2=-Q2/(2D0*m22)*(1D0-BETA)
        C0_Wdec=-1D0/Q2/BETA*(XSPENZ_Wdec(1D0/Y1)-XSPENZ_Wdec(1D0/Y2)
     &                   +1D0/2*(LOG(-Y2)**2-LOG(-Y1)**2)
     &                   +2D0*(LOG(-Y1)-LOG(Y1))*LOG(Y2))
        RETURN
       ELSEIF(rmb2.GE.1D-2.AND.rm22.LT.1D-2) THEN
* ZG light mass case
        ALS=LOG(-Q2/m22)
        C0_Wdec=1D0/Q2*(ALS*LOG((Q2+mb2)/mb2)-zet2
     &         +XSPENZ_Wdec((mb2+Q2)/mb2))
        RETURN
       ELSEIF(rmb2.GE.1D-2.AND.rm22.GE.1D-2) THEN
* ZG massive case is computed below
        GO TO 100
       ELSE
        PRINT *,'Something wrong in c0_sp'
        STOP
       ENDIF
      ENDIF
*
 100  CONTINUE
*
* Cases p12=p22=0
*
      IF(p12.EQ.0D0.AND.p22.EQ.0D0.AND.rm12.EQ.rm32) THEN     
*
* Case p12=p22=0 and m12,m22,m12
*                          
       SQR=SQRT(1D0+4D0*m12/Q2)
       X1=(1D0-SQR)/2D0
       X2=(1D0+SQR)/2D0
       X0=(m22-m12)/Q2
       X3=m22/(m22-m12)
       Y1=X1/X0
       Y2=(1D0-X1)/(1D0-X0)
       Y3=X2/X0
       Y4=(1D0-X2)/(1D0-X0)
       Y5=X3/X0
       Y6=(1D0-X3)/(1D0-X0)
*
       C0_Wdec=1D0/Q2*(XSPENZ_Wdec(1D0/(1D0-Y1))
     &           -XSPENZ_Wdec(1D0/(1D0-Y2))
     &           +XSPENZ_Wdec(1D0/(1D0-Y3))-XSPENZ_Wdec(1D0/(1D0-Y4))
     &           -XSPENZ_Wdec(1D0/(1D0-Y5))+XSPENZ_Wdec(1D0/(1D0-Y6)))
       RETURN
*
      ELSEIF(p12.EQ.0D0.AND.p22.EQ.0D0) THEN     
*
* Case p12=p22=0 and arbitrary m12,m22,m32
*
       LAMBDA=Q2**2+2D0*Q2*(m12+m32)+(m12-m32)**2
       SQR=SQRT(LAMBDA)
       SQS=Q2+m12-m32
       X0=1D0+(m12-m22)/Q2
       X1=(SQS-SQR)/(2D0*Q2)      
       X2=(SQS+SQR)/(2D0*Q2)
       X3=m32/(m32-m22)
       C0_Wdec=1D0/Q2*(XSPENZ_Wdec((X0-1D0)/(X0-X1))
     &           -XSPENZ_Wdec(X0/(X0-X1))
     &           +XSPENZ_Wdec((X0-1D0)/(X0-X2))-XSPENZ_Wdec(X0/(X0-X2))
     &           -XSPENZ_Wdec((X0-1D0)/(X0-X3))+XSPENZ_Wdec(X0/(X0-X3)))
       RETURN
      ENDIF
*
* Pinch of WW-direct box, cat0wb = C0WWB(0d0,-rtm2,-t,cz2,wm2,bm2)
*                         cat0wb = C0WWB(-rtm2,0d0,-t,bm2,wm2,cz2)
*                         cat0w0 = C0WWB(0d0,-rtm2,-t,cz2,wm2,cz2)
*
       IF((p12.EQ.0D0.AND.-p22.NE.rm32.AND.Q2.NE.0D0.AND.rm12.EQ.0D0)
     &.OR.(p22.EQ.0D0.AND.-p12.NE.rm12.AND.Q2.NE.0D0.AND.rm32.EQ.0D0)
     &   ) THEN
       IF(p12.EQ.0D0.AND.rm12.EQ.0D0) THEN
        TMI=p22-Q2
        YPL=(rm22-p22-rm32+DSQRT((rm22-p22-rm32)**2+4D0*rm22*p22))
     &                                                 /2D0/(-p22)
        YMI=rm22/(-p22)/YPL
        YDE=1D0+rm32/Q2
       ENDIF
       IF(p22.EQ.0D0.AND.rm32.EQ.0D0) THEN
        TMI=p12-Q2
        YPL=(rm22-p12-rm12+DSQRT((rm22-p12-rm12)**2+4D0*rm22*p12))
     &                                                 /2D0/(-p12)
        YMI=rm22/(-p12)/YPL
        YDE=1D0+rm12/Q2
       ENDIF
       C0_Wdec=DCMPLX(1D0/TMI*(-DDILOG_Wdec((rm22+TMI)/(rm22+TMI*YPL))
     &                    +DDILOG_Wdec((rm22    )/(rm22+TMI*YPL))
     &                    -DDILOG_Wdec((rm22+TMI)/(rm22+TMI*YMI))
     &                    +DDILOG_Wdec((rm22    )/(rm22+TMI*YMI))
     &                    +DDILOG_Wdec((rm22+TMI)/(rm22+TMI*YDE))
     &                    -DDILOG_Wdec((rm22    )/(rm22+TMI*YDE))
     &                    +DDILOG_Wdec(1D0+TMI/rm22)-zet2),0D0  )
       RETURN
*
* Two W-decay cases
* 
      ELSEIF(-p12.EQ.rm22.AND.p22.EQ.0D0.AND.-Q2.EQ.rm32) THEN
*
       R21=m22/m12
       R31=m32/m12
       BETA21=SQRT(1D0-4D0*R21)
       BETA31=SQRT(1D0-4D0*R31)
       C0_Wdec=1D0/(m32-m22)*(LOG((1D0+BETA31)/2D0)*LOG(2D0*R31
     &         /(1D0+BETA31))
     &                 -LOG((1D0+BETA21)/2D0)*LOG(2D0*R21/(1D0+BETA21)))
       RETURN
* 
      ELSEIF(p12.EQ.0D0.AND.-p22.EQ.rm22.AND.-Q2.EQ.rm12) THEN
*
       R13=m12/m32
       R23=m22/m32
       BETA13=SQRT(1D0-4D0*R13)
       BETA23=SQRT(1D0-4D0*R23)
       C0_Wdec=1D0/(m12-m22)*(LOG((1D0+BETA13)/2D0)*LOG(2D0*R13
     &        /(1D0+BETA13))
     &                 -LOG((1D0+BETA23)/2D0)*LOG(2D0*R23/(1D0+BETA23)))
       RETURN
      ENDIF
*
* Variations of exact case
*
* 1) c0(0d0,-rtm2,-t,cz2,zm2,tm2)
*
       IF((p12.EQ.0D0.AND.-p22.EQ.rm32.AND.Q2.NE.0D0.AND.rm12.EQ.0D0)
     &.OR.(p12.EQ.0D0.AND.-p22.EQ.rm12.AND.Q2.NE.0D0.AND.rm32.EQ.0D0)
     &   ) THEN
*
* Case p12.eq.0.and.(-p22.eq.rm12.or.rm32).and.Q2.ne.0
*
       BETAZ=SQRT(1D0+4D0*p22/m22)
       Y0=m22/(Q2-p22)
       Y1=2D0/(1D0-BETAZ)
       Y2=2D0/(1D0+BETAZ)
       C0_Wdec=1D0/(Q2-p22)*(zet2-XSPENZ_Wdec(1D0-1D0/Y0)
     &                     -XSPENZ_Wdec((1D0-Y0)*Q2/(-p22+Q2-Y0*Q2))
     &                     +XSPENZ_Wdec((   -Y0)*Q2/(-p22+Q2-Y0*Q2))
     &  +XSPENZ_Wdec((1D0-Y0)/(Y1-Y0))-XSPENZ_Wdec((   -Y0)/(Y1-Y0))
     &  +XSPENZ_Wdec((1D0-Y0)/(Y2-Y0))-XSPENZ_Wdec((   -Y0)/(Y2-Y0)))
       RETURN
*
* 2) c0(-rtm2,0d0,-t,tm2,zm2,cz2)
*
      ELSEIF((-p12.EQ.rm12.AND.p22.EQ.0D0.AND.Q2.NE.0D0.AND.rm32.EQ.0D0)
     &   .OR.(-p12.EQ.rm32.AND.p22.EQ.0D0.AND.Q2.NE.0D0.AND.rm12.EQ.0D0)
     &      ) THEN
*
* Case (-p12.eq.rm12.or.rm32).and.p22.eq.0.and.Q2.ne.0
*
       BETAZ=SQRT(1D0+4D0*p12/m22)
       Y0=m22/(Q2-p12)
       Y1=2D0/(1D0-BETAZ)
       Y2=2D0/(1D0+BETAZ)
       C0_Wdec=1D0/(Q2-p12)*(zet2-XSPENZ_Wdec(1D0-1D0/Y0)
     &                     -XSPENZ_Wdec((1D0-Y0)*Q2/(-p12+Q2-Y0*Q2))
     &                     +XSPENZ_Wdec((   -Y0)*Q2/(-p12+Q2-Y0*Q2))
     &  +XSPENZ_Wdec((1D0-Y0)/(Y1-Y0))-XSPENZ_Wdec((   -Y0)/(Y1-Y0))
     &  +XSPENZ_Wdec((1D0-Y0)/(Y2-Y0))-XSPENZ_Wdec((   -Y0)/(Y2-Y0)))
       RETURN
* 3)
      ELSEIF(p12.NE.0D0.AND.p22.EQ.0D0.AND.-Q2.EQ.rm12.AND.rm22.EQ.0D0
     &      ) THEN
*
* Case p12.ne.0.and.-Q2.eq.rm12.and.p2.eq.0
*
       BETAZ=SQRT(1D0+4D0*Q2/m32)
       Y0=m32/(p12-Q2)
       Y1=2D0/(1D0-BETAZ)
       Y2=2D0/(1D0+BETAZ)
       C0_Wdec=1D0/(p12-Q2)*(zet2-XSPENZ_Wdec(1D0-1D0/Y0)
     &                     -XSPENZ_Wdec((1D0-Y0)*p12/(-Q2+p12-Y0*p12))
     &                     +XSPENZ_Wdec((   -Y0)*p12/(-Q2+p12-Y0*p12))
     &  +XSPENZ_Wdec((1D0-Y0)/(Y1-Y0))-XSPENZ_Wdec((   -Y0)/(Y1-Y0))
     &  +XSPENZ_Wdec((1D0-Y0)/(Y2-Y0))-XSPENZ_Wdec((   -Y0)/(Y2-Y0)))
       RETURN
* 4)
      ELSEIF(p12.EQ.0D0.AND.p22.NE.0D0.AND.-Q2.EQ.rm32.AND.rm22.EQ.0D0
     &      ) THEN
*
* Case -p12.eq.0.and.p22.ne.0.and.Q2.eq.rm32; PROBLEMS HERE
*
       BETAZ=SQRT(1D0+4D0*Q2/m12)
       Y0=m12/(p22-Q2)
       Y1=2D0/(1D0-BETAZ)
       Y2=2D0/(1D0+BETAZ)
       C0_Wdec=1D0/(p22-Q2)*(zet2-XSPENZ_Wdec(1D0-1D0/Y0)
     &                     -XSPENZ_Wdec((1D0-Y0)*p22/(-Q2+p22-Y0*p22))
     &                     +XSPENZ_Wdec((   -Y0)*p22/(-Q2+p22-Y0*p22))
     &  +XSPENZ_Wdec((1D0-Y0)/(Y1-Y0))-XSPENZ_Wdec((   -Y0)/(Y1-Y0))
     &  +XSPENZ_Wdec((1D0-Y0)/(Y2-Y0))-XSPENZ_Wdec((   -Y0)/(Y2-Y0)))
       RETURN
* 5)
      ELSEIF(p12.NE.0D0.AND.p22.NE.0D0) THEN
*
* Case p12.ne.0.and.p22.ne.0; or GENERAL case  
*
      IMI =DCMPLX(0D0,1D0)
      p1dp2=(Q2-p12-p22)/2D0
      DET3 =p1dp2**2-p12*p22
      SQD3 =DSQRT(DET3)
      IF(p1dp2.ge.0d0) THEN
       alpha=(p1dp2+SQD3)/(-p12)
      ELSE
       alpha=(-p22)/(p1dp2-SQD3)
      ENDIF
      N=2D0*SQD3
c     print *,'alpha,N=',alpha,N
      d=DREAL(p22+m22-m32)
      e=DREAL(Q2-p22+m12-m22)
      r0=(-d-e*alpha)/N
*
      YR(0) =-r0/alpha
      YR(1) = r0/(1d0-alpha)
      YR(2) = r0+alpha
c     print *,'Roots=',YR(0),YR(1),YR(2)
*
      AP(0) =-p22
      BP(0) = p22+m22-m32
      CP(0) = m32
*       
      AP(1) =-Q2
      BP(1) = Q2+m12-m32
      CP(1) = m32
*  
      AP(2) =-p12
      BP(2) = p12+m12-m22
      CP(2) = m22
*
      EP(0)=eps
      EP(1)=eps
      EP(2)=eps
*
      DE(0)=-DIMAG(AP(0)*YR(0)**2+BP(0)*YR(0)+CP(0)-IMI*eps)
      DE(1)=-DIMAG(AP(1)*YR(1)**2+BP(1)*YR(1)+CP(1)-IMI*eps)
      DE(2)=-DIMAG(AP(2)*YR(2)**2+BP(2)*YR(2)+CP(2)-IMI*eps)
*
      C0_Wdec=0D0
*
      DO I=0,2,1
      DP(I)=BP(I)**2-4D0*AP(I)*CP(I)
c     print *,'I DP(I)=',I,DP(I)
      SQDP(I)=SQRT(DP(I))      
      YP(I)=(-BP(I)+SQDP(I))/2D0/AP(I)
c     print *,'I YP(I)=',I,YP(I)
      IF(DABS(DREAL(YP(I))).LT.1D-7) THEN
       YP(I)=-2D0*CP(I)/(BP(I)+SQDP(I))
c      print *,'I YP(I)=',I,YP(I)
      ENDIF
      YM(I)=(-BP(I)-SQDP(I))/2D0/AP(I)
c     print *,'I YM(I)=',I,YM(I)
      IF(DABS(DREAL(YM(I))).LT.1D-7) THEN
       YM(I)=-2D0*CP(I)/(BP(I)-SQDP(I))
c      print *,'I YM(I)=',I,YM(I)
      ENDIF
      IF(Q2.GT.0D0) THEN
        RLOG=DLOG(DABS(1D0-1D0/YR(I))) 
      ELSE
        RLOG=DLOG(DABS(1D0-1D0/YR(I))) 
      ENDIF
      SKERN1(I)=-(+VELTET_Wdec(-YP(I),-YM(I))
     &         -VELTET_Wdec(YR(I)-YP(I),YR(I)-YM(I))
     &           )*RLOG
      SKERN2(I)=-(-VELTET_Wdec(AP(I)-IMI*EP(I),1D0/(AP(I)-IMI*DE(I))))
     &         *RLOG
      SKERN(I)=SKERN1(I)+SKERN2(I)
*
* Analyser SKERN1
*
*     IF(DREAL(SKERN1(I)).NE.0D0.OR.DIMAG(SKERN1(I)).NE.0D0) 
*    &              PRINT *,'I,Re(SKERN1(I)),Im(SKERN1(I))=',I,SKERN1(I)
*
* Analyser SKERN2
*
*     IF(DREAL(SKERN2(I)).NE.0D0.OR.DIMAG(SKERN2(I)).NE.0D0) 
*    &              PRINT *,'I,Re(SKERN2(I)),Im(SKERN2(I))=',I,SKERN2(I)
*
      TERM=(-1)**I/N*(RKERN_Wdec(YR(I),YP(I))+RKERN_Wdec(YR(I),YM(I))
     &    +SKERN(I))
      C0_Wdec=C0_Wdec+TERM
      ENDDO
*
      ELSE
       PRINT *,'Combination of arguments is not yet foreseen!'
       STOP
      ENDIF
*
      RETURN
      END

      COMPLEX*16 FUNCTION RKERN_Wdec(Y,YL)
*
      IMPLICIT NONE
      REAL*8 Y,DREAL,DIMAG
      COMPLEX*16 YL,RKERNR,RKERNI,XSPENZ_Wdec,LOG,VELTET_Wdec
*
      RKERNR=XSPENZ_Wdec(Y/(Y-YL))-XSPENZ_Wdec((Y-1D0)/(Y-YL))
      RKERNI=VELTET_Wdec(   -YL,1D0/(Y-YL))*LOG(Y/(Y-YL))
     &      -VELTET_Wdec(1D0-YL,1D0/(Y-YL))*LOG((Y-1D0)/(Y-YL))
      RKERN_Wdec=RKERNR+RKERNI
*
* Analyser RKERN
*
      IF(DREAL(RKERNI).NE.0D0.OR.DIMAG(RKERNI).NE.0D0) 
     &                                   PRINT *,'Re(RI),Im(RI)=',RKERNI
*
      RETURN
      END

      COMPLEX*16 FUNCTION VELTET_Wdec(A,B)
*
      IMPLICIT NONE
      REAL*8 PI,IMA,IMB,IMAB,TET_Wdec

      COMPLEX*16 A,B,IMI,DCMPLX
*
      DATA PI/3.14159274101D0/
*
      IMA =DIMAG(A)
      IMB =DIMAG(B)
      IMAB=DIMAG(A*B)
      IMI =DCMPLX(0D0,1D0)
      VELTET_Wdec=2D0*IMI*PI
     &           *(TET_Wdec(-IMA)*TET_Wdec(-IMB)*TET_Wdec( IMAB)
     &                  -TET_Wdec( IMA)*TET_Wdec( IMB)*TET_Wdec(-IMAB))
*
      RETURN
      END

      REAL*8 FUNCTION TET_Wdec(X)
*
      IMPLICIT NONE
      REAL*8 X
*
      IF(X.LT.0D0) THEN
        TET_Wdec=0D0
      ELSE
        TET_Wdec=1D0
      ENDIF
*
      RETURN
      END
