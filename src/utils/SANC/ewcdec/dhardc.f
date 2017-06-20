      real*8 function dhardp2_Wdec(ssp)
      IMPLICIT NONE
      real*8 ss,cc,sp,v2,em2,fm2,reps,aeps,ssp,dhardp_Wdec
      common/hard/ss,cc,sp,v2,em2,fm2,reps,aeps
      dhardp2_Wdec= dhardp_Wdec(ss,ssp)
      return
      end

c begin of package hard ============================================
      real*8 function dhardc_Wdec(s,cost,omega)
c ====================================================
c package hard Arnd Leike 19.12.00
c calculates O(alpha**3) hard photon corrections to e+e- --> tt(bar)
c input:  s = c.m. energy squared in GeV**2, 
c         cost = cos( theta (e+,t(bar) ) ) 
c         omega = minimum photon energy in GeV
c output: d sigma / d cost  in  pb
c
c  numerical integration of fsp over sp
c
c Here the first version with 2 numerical integrations over v2 and sp.
c After succesfull tests, the numerical integrations have to be replaced
c by (partly existing) analytical expressions.
c
c contains the following subroutines:
c dhards: integrates fsp over sp gives the output described above
c fsp: integrates fv2ini, fv2int, fv2fin over v2
c fv2ini,fv2int,fv2fin: squared matrix element integrated over phig
c                       analytically for ini, int, fin
c integration scheme of the first version:
c fv2ini,fv2int,fv2fin ==> over_v2 ==> fsp ==> over_sp ==> dhardc
      IMPLICIT NONE
*-
      real*8 s,cost,omega,delta,adelta,a,b,h,out,aih,aiabs,fsp_Wdec
      REAL*8 rmf,Qf,cf,I3f
      real*8      ss,cc,sp,v2,em2,fm2,reps,aeps
      common/hard/ss,cc,sp,v2,em2,fm2,reps,aeps
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
c       nue,e,numu,mu,nutau,tau,um,dm,cm,sm,tm,bm
      external fsp_Wdec
c     fill common hard as far as possible
      ss=s
      cc=cost
      em2=rmf(2)**2
      fm2=rmf(11)**2
      reps=1d-4
      aeps=0d0
      a=4d0*fm2*(1d0+1d-6)
c     addition of 1d-6 for numerical reasons
c delta: cut in the photon energy changes the lower integration bound a
      delta=0.5d0
      adelta=ss*(1d0-delta) 
c      if(a.lt.adelta) a=adelta
      b=ss-2d0*sqrt(s)*omega
      h=(b-a)/3d0
      write(6,100) a,b,h,em2,fm2
 100  format(' hardin:',5e12.4) 
      call SIMPS1_Wdec(a,b,h,reps,aeps,fsp_Wdec,sp,out,aih,aiabs)
      dhardc_Wdec=out
      return
      end

      real*8 function fsp_Wdec(ssp)
c ====================================================
c fsp is part of the package hard:  Arnd Leike 15.12.00
c input: ssp = sp = ( p(top) + p(top(bar)) )^2
c output: d^2 sigma / d cost d sp  in  pb/GeV
c
c  numerical integration of fv2ini, fv2int or fv2fin over v2
c
      IMPLICIT NONE
      real*8 ssp,a,b,h,out,aih,aiabs,
     &       fv2int_Wdec,fv2ini_Wdec,fv2fin_Wdec_Wdec,dsqrt
      real*8      ss,cc,sp,v2,em2,fm2,reps,aeps
      common/hard/ss,cc,sp,v2,em2,fm2,reps,aeps
      external fv2int_Wdec,fv2ini_Wdec,fv2fin_Wdec
      sp=ssp
      a=(ss-sp)*(1d0-dsqrt(1d0-4d0*fm2/sp))/2d0*(1d0+1d-6)
      b=(ss-sp)*(1d0+dsqrt(1d0-4d0*fm2/sp))/2d0*(1d0-1d-6)
c     additions of +-1d-6 for numerical reasons 
      h=(b-a)/3d0
c      write(6,100) a,b,h,sp
 100  format(' fspin:',5e12.4) 
      call SIMPS2_Wdec(a,b,h,reps,aeps,fv2ini_Wdec,v2,out,aih,aiabs)
c      call SIMPS2_Wdec(a,b,h,reps,aeps,fv2int_Wdec,v2,out,aih,aiabs)
c      call SIMPS2_Wdec(a,b,h,reps,aeps,fv2fin_Wdec,v2,out,aih,aiabs)
      fsp_Wdec=out
c
c      write(6,1100) cc,sp,dsqrt(sp),fsp
 1100 format(' fsp:',4e16.8)
      return
      end

      real*8 function dhardp_Wdec(s,ssp)
c ====================================================
c dhardp is part of the package hard Arnd Leike 19.12.00
c calculates O(alpha**3) hard photon corrections to e+e- --> tt(bar)
c input:  s = c.m. energy squared in GeV**2, 
c         sp = ( p(t) + p(t(bar)) )**2 
c output: d sigma / d sp  in  pb/GeV
c
c  numerical integration of fsp over cost
c integration scheme 
c fv2ini,fv2int,fv2fin ==> over_v2 ==> fcost ==> over_cost ==> dhardp
      IMPLICIT NONE
*-
      real*8 s,ssp,cost,a,b,h,out,aih,aiabs,fcost_Wdec
      REAL*8 rmf,Qf,cf,I3f
      real*8      ss,cc,sp,v2,em2,fm2,reps,aeps
      common/hard/ss,cc,sp,v2,em2,fm2,reps,aeps
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
c       nue,e,numu,mu,nutau,tau,um,dm,cm,sm,tm,bm
      data a,b,h/-0.999999d0,0.999999d0,0.6d0/
      external fcost_Wdec
c     fill common hard as far as possible
      ss=s
      sp=ssp
      em2=rmf(2)**2
      fm2=rmf(11)**2
c      fm2=1d0**2
      reps=1d-4
      aeps=0d0
      call SIMPS1_Wdec(a,b,h,reps,aeps,fcost_Wdec,cost,out,aih,aiabs)
      dhardp_Wdec=out
      write(6,200) sqrt(s),sqrt(sp),out*2d0*sqrt(sp)
 200     format(' sqrt(s),sqrt(sp),out', 4d12.4)
      return
      end

      real*8 function fcost_Wdec(cost)
c ====================================================
c fcost is part of the package hard:  Arnd Leike 20.12.00
c input: cost = cos( theta (e+,t(bar) ) ) 
c output: d^2 sigma / d cost d sp  in  pb/GeV
c
c  numerical integration of fv2ini, fv2int or fv2fin over v2
c
      IMPLICIT NONE
      real*8 cost,a,b,h,out,aih,aiabs
      real*8 fv2int_Wdec,fv2ini_Wdec,fv2fin_Wdec,dsqrt
      real*8      ss,cc,sp,v2,em2,fm2,reps,aeps
      common/hard/ss,cc,sp,v2,em2,fm2,reps,aeps
      external fv2int_Wdec,fv2ini_Wdec,fv2fin_Wdec
      cc=cost
      a=(ss-sp)*(1d0-dsqrt(1d0-4d0*fm2/sp))/2d0*(1d0+1d-6)
      b=(ss-sp)*(1d0+dsqrt(1d0-4d0*fm2/sp))/2d0*(1d0-1d-6)
c     additions of +-1d-6 for numerical reasons 
      h=(b-a)/3d0
c      write(6,100) a,b,h,sp
 100  format(' fspin:',5e12.4) 
      call SIMPS2_Wdec(a,b,h,reps,aeps,fv2ini_Wdec,v2,out,aih,aiabs)
c      call SIMPS2_Wdec(a,b,h,reps,aeps,fv2int_Wdec,v2,out,aih,aiabs)
c      call SIMPS2_Wdec(a,b,h,reps,aeps,fv2fin_Wdec,v2,out,aih,aiabs)
      fcost_Wdec=out
c
c      write(6,1100) cc,sp,dsqrt(sp),fsp
 1100 format(' fsp:',4e16.8)
      return
      end

      real*8 function fv2ini_Wdec(vv2)
c ====================================================
c fv2ini is part of the package hard:  Arnd Leike 19.12.00
c input: vv2 = v2 = -2 * p(gamma) * p(top(bar))
c output: d^3 sigma / d cost d sp  d v2  in  pb/GeV^2
c
c  analytical integration of the squared matrix element over phig
c initial state radiation
c
      IMPLICIT NONE
      INTEGER*4 ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      real*8 v1,tt,uu,ctg,beta0,dsga,dspga,pi,kappa,const3
      real*8 vv2,slam,c1,sqc1,c2,sqc2,z1,z2,z1i,z2i,z1i2,z2i2,z12i
      real*8 vfactp,afactp,cfactp,af,vf,qff,gamma
      real*8 dsqrt,dreal
      real*8      ss,cc,sp,v2,em2,fm2,reps,aeps
      REAL*8 rmf,Qf,cf,I3f
      real*8 CONHC,alphaI,GFermi,GammaZ,DAL5H
      real*8 zm,rwm2,rzm2,rhm2
      real*8 qt,qtm,qb,qe,anu,vnu,ae,ve,ae2,ve2,AAe,qe2,qt2,qeqt,
     &          at,vt,at2,vt2,AAt,AAv,ab,vb,vpae,vmae,vpau,vmau,de,dt,Nc
      COMPLEX*16 xmz2,xsz,xszg,xspz,xspzg,cz2,wm2,zm2,hm2,xkapp
      complex*16 dcmplx,dconjg
c
      common/hard/ss,cc,sp,v2,em2,fm2,reps,aeps
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
c       nue,e,numu,mu,nutau,tau,um,dm,cm,sm,tm,bm
      COMMON/FLAGGS/ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      COMMON/PHYSCO/CONHC,alphaI,GFermi,GammaZ,DAL5H
      COMMON/BOS_m2/zm,cz2,wm2,zm2,hm2,rwm2,rzm2,rhm2
      COMMON/COUPLC/qt,qtm,qb,qe,anu,vnu,ae,ve,ae2,ve2,AAe,qe2,qt2,qeqt,
     &          at,vt,at2,vt2,AAt,AAv,ab,vb,vpae,vmae,vpau,vmau,de,dt,Nc
c common calculations:
      v2=vv2
      PI = 3.141592654D0
      const3=conhc/alphaI**3/4d0/ss**2/pi
      xmz2=dcmplx(rzm2,-zm*GammaZ)
      kappa=GFermi*rzm2/(dSQRT(2d0)*2d0*PI)*alphaI
c This kappa is 4 times larger than in ZFITTER because the Normalization
c of the Z-couplings to fermions is different in TOPFIT and ZFITTER
      IF(GAMS.EQ.0) THEN
        xkapp =kappa
          ELSE
c Z-Boson parameter transformation
        gamma =GammaZ/zm
        xmz2  =xmz2 /(1d0+gamma**2)
        xkapp =kappa/(1d0+(0d0,1d0)*gamma)
      ENDIF
c normalized propagators
       dspga=1d0
       xspz=(sp-xmz2)/xkapp/sp
       xspzg=dconjg(xspz)
       dsga=1d0
       xsz=(ss-xmz2)/xkapp/ss
       xszg=dconjg(xsz)
c some invariants
      v1=ss-sp-v2
      tt=(ss-v1 + cc*dsqrt((ss-v1)**2-4d0*fm2*ss))/2d0
      uu=sp+v2-tt
c table of integrals
      ctg=(v1*(ss+sp)-ss*(ss-sp))/(ss-sp)/dsqrt((ss-v1)**2-4d0*fm2*ss)
      beta0=dsqrt(1d0-4d0*em2/ss)
c
      slam=dsqrt((sp+v2)**2-4d0*fm2*ss)
      c1=(2d0*ss*sp-(v2+sp)*(ss+sp)+(ss-sp)*slam*beta0*cc)**2/4d0
     &   +4d0*em2*(sp*v2*(ss-sp-v2)-(ss-sp)**2*fm2)
      sqc1=dsqrt(c1)
      c2=(2d0*ss*sp-(v2+sp)*(ss+sp)-(ss-sp)*slam*beta0*cc)**2/4d0
     &   +4d0*em2*(sp*v2*(ss-sp-v2)-(ss-sp)**2*fm2)
      sqc2=dsqrt(c2)
      z1=(ss-sp)*(1d0+beta0*cc*ctg)/2d0
      z2=(ss-sp)*(1d0-beta0*cc*ctg)/2d0
      z1i=slam/sqc1
      z2i=slam/sqc2
      z1i2=z1*z1i**3
      z2i2=z2*z2i**3
      z12i=(z1i+z2i)/(z1+z2)
c copulings
      qff=qt
      vf=vt
      af=at
c
      vfactp= qe**2*qff**2/dspga**2 
     &      + 2d0*ve*qe*vf*qff*dreal(1/xspz/dspga)
     &      + (ve**2+ae**2)*(vf**2+af**2)*dreal(1/xspz/xspzg) 
      afactp= 
     &      + 2d0*ae*qe*af*qff*dreal(1/xspz/dspga)
     &      + (ve*ae+ae*ve)*(vf*af+af*vf)*dreal(1/xspz/xspzg) 
      cfactp=4d0*(ve**2+ae**2)*af**2*dreal(1/xspz/xspzg)  
c
      fv2ini_Wdec=const3/sp**2*2d0*pi*(
     & Vfactp*(-2*em2*Z1i2*(2*UU**2-2*UU*SP+SP**2+2*fm2*SP)
     &         -2*em2*Z2i2*(2*TT**2-2*TT*SP+SP**2+2*fm2*SP)
     &+SP*Z12i*(2*UU**2-2*UU*SP+2*TT**2-2*TT*SP+2*SP**2+4*fm2*SP)
     &+Z1i*(-2*TT*SP+SS*SP+SP**2+2*fm2*(SS+SP))
     &+Z2i*(-2*UU*SP+SS*SP+SP**2+2*fm2*(SS+SP))  -2*SP-4*fm2 )
     &+     Afactp*SP*(-2*em2*Z1i2*(-2*UU+SP) +2*em2*Z2i2*(-2*TT+SP)
     &+2*SP*Z12i*(TT-UU) -Z1i*(-2*TT+SS+SP) +Z2i*(-2*UU+SS+SP)  )
     &+     Cfactp*fm2*(2*em2*Z1i2*SP +2*em2*Z2i2*SP 
     &-2*SP**2*Z12i - (SS+SP)*(Z1i+Z2i) +2 )   )           
      return
      end

      real*8 function fv2int_Wdec(vv2)
c ====================================================
c fv2int is part of the package dhardc:  Arnd Leike 19.12.00
c input: vv2 = v2 = -2 * p(gamma) * p(top(bar))
c output: d^3 sigma / d cost d sp  d v2  in  pb/GeV^2
c
c  analytical integration of the squared matrix element over phig
c interference between initial and final state radiation
      IMPLICIT NONE
      INTEGER*4 ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      real*8 v1,tt,uu,ctg,beta0,z1,z2,dsga,dspga,pi,kappa,const3
      real*8 vv2,slam,c1,sqc1,c2,sqc2,z1i,z2i
      real*8 vfap,afap,cfap,cfapmi,af,vf,qff,gamma
      real*8 dsqrt,dreal
      real*8      ss,cc,sp,v2,em2,fm2,reps,aeps
      REAL*8 rmf,Qf,cf,I3f
      real*8 CONHC,alphaI,GFermi,GammaZ,DAL5H
      real*8 zm,rwm2,rzm2,rhm2
      real*8 qt,qtm,qb,qe,anu,vnu,ae,ve,ae2,ve2,AAe,qe2,qt2,qeqt,
     &          at,vt,at2,vt2,AAt,AAv,ab,vb,vpae,vmae,vpau,vmau,de,dt,Nc
      COMPLEX*16 xmz2,xsz,xszg,xspz,xspzg,cz2,wm2,zm2,hm2,xkapp
      complex*16 dcmplx,dconjg
c
      common/hard/ss,cc,sp,v2,em2,fm2,reps,aeps
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
c       nue,e,numu,mu,nutau,tau,um,dm,cm,sm,tm,bm
      COMMON/FLAGGS/ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      COMMON/PHYSCO/CONHC,alphaI,GFermi,GammaZ,DAL5H
      COMMON/BOS_m2/zm,cz2,wm2,zm2,hm2,rwm2,rzm2,rhm2
      COMMON/COUPLC/qt,qtm,qb,qe,anu,vnu,ae,ve,ae2,ve2,AAe,qe2,qt2,qeqt,
     &          at,vt,at2,vt2,AAt,AAv,ab,vb,vpae,vmae,vpau,vmau,de,dt,Nc
c common calculations:
      v2=vv2
      PI = 3.141592654D0
      const3=conhc/alphaI**3/4d0/ss**2/pi
      xmz2=dcmplx(rzm2,-zm*GammaZ)
      kappa=GFermi*rzm2/(dSQRT(2d0)*2d0*PI)*alphaI
c This kappa is 4 times larger than in ZFITTER because the Normalization
c of the Z-couplings to fermions is different in TOPFIT and ZFITTER
      IF(GAMS.EQ.0) THEN
        xkapp =kappa
          ELSE
c Z-Boson parameter transformation
        gamma =GammaZ/zm
        xmz2  =xmz2 /(1d0+gamma**2)
        xkapp =kappa/(1d0+(0d0,1d0)*gamma)
      ENDIF
c normalized propagators
       dspga=1d0
       xspz=(sp-xmz2)/xkapp/sp
       xspzg=dconjg(xspz)
       dsga=1d0
       xsz=(ss-xmz2)/xkapp/ss
       xszg=dconjg(xsz)
c some invariants
      v1=ss-sp-v2
      tt=(ss-v1 + cc*dsqrt((ss-v1)**2-4d0*fm2*ss))/2d0
      uu=sp+v2-tt
      ctg=(v1*(ss+sp)-ss*(ss-sp))/(ss-sp)/dsqrt((ss-v1)**2-4d0*fm2*ss)
      beta0=dsqrt(1d0-4d0*em2/ss)
c
      slam=dsqrt((sp+v2)**2-4d0*fm2*ss)
      c1=(2d0*ss*sp-(v2+sp)*(ss+sp)+(ss-sp)*slam*beta0*cc)**2/4d0
     &   +4d0*em2*(sp*v2*(ss-sp-v2)-(ss-sp)**2*fm2)
      sqc1=dsqrt(c1)
      c2=(2d0*ss*sp-(v2+sp)*(ss+sp)-(ss-sp)*slam*beta0*cc)**2/4d0
     &   +4d0*em2*(sp*v2*(ss-sp-v2)-(ss-sp)**2*fm2)
      sqc2=dsqrt(c2)
      z1=(ss-sp)*(1d0+beta0*cc*ctg)/2d0
      z2=(ss-sp)*(1d0-beta0*cc*ctg)/2d0
      z1i=slam/sqc1
      z2i=slam/sqc2
c copulings
      qff=qt
      vf=vt
      af=at
c
      vfap= qe**2*qff**2/dspga/dsga
     &+ ve*qe*vf*qff*dreal(1d0/dspga/xsz+1d0/dsga/xspz)
     &+(ve**2+ae**2)*(vf**2+af**2)*dreal(1d0/xsz/xspzg)
      afap= 
     &  qe*qff*ae*af*dreal(1d0/dspga/xsz+1d0/dsga/xspz)
     &+4d0*ve*ae*vf*af*dreal(1d0/xsz/xspzg)
      cfap=4d0*(ve**2+ae**2)*af**2*dreal(1/xspz/xszg) 
      cfapmi= 2d0*ae*qe*af*qff*dreal(1d0/dspga/xsz-1d0/dsga/xspz)
c
      fv2int_Wdec=const3/ss/sp*2d0*pi*(
     & vfap*((SS-TT)*Z1i/V1*(2*UU**2-2*UU*SP+SP**2+2*fm2*SP
     &-2*TT*UU+SS**2+2*fm2*SS) -(SS-UU)*Z2i/V1*(2*TT**2-2*TT*SP+SP**2
     &+2*fm2*SP -2*TT*UU+SS**2+2*fm2*SS) +(SP-TT)*Z2i/V2*(2*UU**2
     &-2*UU*SS+SS**2+2*fm2*SS -2*TT*UU+SP**2+2*fm2*SP) -(SP-UU)*Z1i/V2*
     &(2*TT**2-2*TT*SS+SS**2+2*fm2*SS -2*TT*UU+SP**2+2*fm2*SP)
     &+1/V1*( (Z1-Z2)*SS          +(TT-UU)*(3*SS-SP-4*fm2) )
     &-1/V2*( (Z1-Z2)*(SP+4*fm2) +(TT-UU)*(3*SP-SS+4*fm2) )
     &       + (Z2i  -Z1i)*(SS**2+SP**2+2*fm2*(SS+SP) )        )
     &+afap*((SS-TT)*Z1i/V1*(-2*UU*SP +SP**2  +SS*(TT-UU)  )
     &        +(SS-UU)*Z2i/V1*(-2*TT*SP +SP**2  -SS*(TT-UU)  )
     &        +(SP-TT)*Z2i/V2*(-2*UU*SS +SS**2  +SP*(TT-UU)  )
     &        +(SP-UU)*Z1i/V2*(-2*TT*SS +SS**2  -SP*(TT-UU)  )
     &+1/V1*(  2*SS*SP+2*fm2*(V2-2*SS) ) +Z1i*(-SS*TT+SP*UU)
     &+1/V2*(- 2*SS*SP+2*fm2*(V1-2*SS) ) +Z2i*(-SS*UU+SP*TT)
     &+2*(SS+SP) + 4*fm2         ) 
     &+cfap*fm2*((SS+SP)*( -(SS-TT)*Z1i/V1 +(SS-UU)*Z2i/V1 
     &                - UU*Z2i/V2 +TT*Z1i/V2 ) +(Z1-Z2)*(1/V2 -1/V1) ) 
     &-cfapmi*fm2*(SS+SP)*(1/V1 +1/V2)  )
      return
      end

      real*8 function fv2fin_Wdec(vv2)
c ====================================================
c fv2fin is part of the package hard:  Arnd Leike 19.12.00
c input: vv2 = v2 = -2 * p(gamma) * p(top(bar))
c output: d^3 sigma / d cost d sp  d v2  in  pb/GeV^2
c
c  analytical integration of the squared matrix element over phig
c final state radiation
      IMPLICIT NONE
      INTEGER*4 ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      real*8 v1,tt,uu,ctg,beta0,dsga,dspga,pi,kappa,const3
      real*8 vv2,slam,c1,sqc1,c2,sqc2,z1,z2,z12
      real*8 vfact,afact,cfact,af,vf,qff,gamma
      real*8 dsqrt,dreal
      real*8      ss,cc,sp,v2,em2,fm2,reps,aeps
      REAL*8 rmf,Qf,cf,I3f
      real*8 CONHC,alphaI,GFermi,GammaZ,DAL5H
      real*8 zm,rwm2,rzm2,rhm2
      real*8 qt,qtm,qb,qe,anu,vnu,ae,ve,ae2,ve2,AAe,qe2,qt2,qeqt,
     &          at,vt,at2,vt2,AAt,AAv,ab,vb,vpae,vmae,vpau,vmau,de,dt,Nc
      COMPLEX*16 xmz2,xsz,xszg,xspz,xspzg,cz2,wm2,zm2,hm2,xkapp
      complex*16 dcmplx,dconjg
c
      common/hard/ss,cc,sp,v2,em2,fm2,reps,aeps
      COMMON/FERMrm/rmf(12),Qf(12),cf(12),I3f(12)
c       nue,e,numu,mu,nutau,tau,um,dm,cm,sm,tm,bm
      COMMON/FLAGGS/ALEM,ALE2,VPOL,QCDC,ITOP,BOXD,DIAG,GAMS,WEAK,IMOMS
      COMMON/PHYSCO/CONHC,alphaI,GFermi,GammaZ,DAL5H
      COMMON/BOS_m2/zm,cz2,wm2,zm2,hm2,rwm2,rzm2,rhm2
      COMMON/COUPLC/qt,qtm,qb,qe,anu,vnu,ae,ve,ae2,ve2,AAe,qe2,qt2,qeqt,
     &          at,vt,at2,vt2,AAt,AAv,ab,vb,vpae,vmae,vpau,vmau,de,dt,Nc
c common calculations:
      v2=vv2
      PI = 3.141592654D0
      const3=conhc/alphaI**3/4d0/ss**2/pi
      xmz2=dcmplx(rzm2,-zm*GammaZ)
      kappa=GFermi*rzm2/(dSQRT(2d0)*2d0*PI)*alphaI
c This kappa is 4 times larger than in ZFITTER because the Normalization
c of the Z-couplings to fermions is different in TOPFIT and ZFITTER
      IF(GAMS.EQ.0) THEN
        xkapp =kappa
          ELSE
c Z-Boson parameter transformation
        gamma =GammaZ/zm
        xmz2  =xmz2 /(1d0+gamma**2)
        xkapp =kappa/(1d0+(0d0,1d0)*gamma)
      ENDIF
c normalized propagators
       dspga=1d0
       xspz=(sp-xmz2)/xkapp/sp
       xspzg=dconjg(xspz)
       dsga=1d0
       xsz=(ss-xmz2)/xkapp/ss
       xszg=dconjg(xsz)
c some invariants
      v1=ss-sp-v2
      tt=(ss-v1 + cc*dsqrt((ss-v1)**2-4d0*fm2*ss))/2d0
      uu=sp+v2-tt
c table of integrals
      ctg=(v1*(ss+sp)-ss*(ss-sp))/(ss-sp)/dsqrt((ss-v1)**2-4d0*fm2*ss)
      beta0=dsqrt(1d0-4d0*em2/ss)
c
      slam=dsqrt((sp+v2)**2-4d0*fm2*ss)
      c1=(2d0*ss*sp-(v2+sp)*(ss+sp)+(ss-sp)*slam*beta0*cc)**2/4d0
     &   +4d0*em2*(sp*v2*(ss-sp-v2)-(ss-sp)**2*fm2)
      sqc1=dsqrt(c1)
      c2=(2d0*ss*sp-(v2+sp)*(ss+sp)-(ss-sp)*slam*beta0*cc)**2/4d0
     &   +4d0*em2*(sp*v2*(ss-sp-v2)-(ss-sp)**2*fm2)
      sqc2=dsqrt(c2)
      z1=(ss-sp)*(1d0+beta0*cc*ctg)/2d0
      z2=(ss-sp)*(1d0-beta0*cc*ctg)/2d0
      z12=z1*z2 - ((ss-sp)*beta0)**2*(1d0-cc*cc)*(1d0-ctg*ctg)/8d0
c copulings
      qff=qt
      vf=vt
      af=at
c
      vfact= qe**2*qff**2/dsga**2 
     &      + 2d0*ve*qe*vf*qff*dreal(1/xsz/dsga)
     &      + (ve**2+ae**2)*(vf**2+af**2)*dreal(1/xsz/xszg) 
      afact= 
     &      + 2d0*ae*qe*af*qff*dreal(1/xsz/dsga)
     &      + (ve*ae+ae*ve)*(vf*af+af*vf)*dreal(1/xsz/xszg) 
      cfact=4d0*(ve**2+ae**2)*af**2*dreal(1/xsz/xszg) 
c
      fv2fin_Wdec=const3/ss/ss*2*pi*(
     &vfact*(-2*fm2/V1**2*(-2*TT*UU +SS**2 +2*fm2*SS)
     &  -2*fm2/V2**2*(-2*(tt*uu+tt*z2+z1*uu+z12)+SS**2 +2*fm2*SS)
     &  +4*fm2/V1/V2*(TT*(UU+Z2)+UU*(TT+Z1)    -2*fm2*SS) 
     &  -SS/V1/V2*( 2*TT*UU +  2*(tt*uu+tt*z2+z1*uu+z12) -2*SS**2) 
     &  +SS/V1*(V2 -4*fm2) +SS/V2*(V1 -2*SS -4*fm2)  )
     &+afact*SS*(2*fm2/V1**2*( -TT+UU ) + 2*fm2/V2**2*(-(TT+Z1)+UU+Z2)
     &+(SP-2*fm2)/V1/V2*(2*TT-2*UU+Z1-Z2)+(TT-UU)/V1+(TT+Z1 -UU-Z2)/V2)
     &+cfact*fm2*(2*(Z12 -SS*SP +2*fm2*SS)/V1/V2 
     &   + 2*fm2*SS*(1/V1**2+1/V2**2)  )   )
      return
      end

      SUBROUTINE SIMPS1_Wdec(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C =====================================================================
C SIMPS
C A1,B1 -THE LIMITS OF INTEGRATION
C H1    -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN  AEPS1=REPS1=0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(7),P(5)
      external funct
      H=DSIGN(H1,B1-A1)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(X)/3.D0
    4 X0=X
      IF((X0+4.D0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.D0
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(X)/3.D0
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*ABS(F(K))
      DI1=(F(1)+4.D0*F(3)+F(5))*2.D0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.D0) 17,14,14
   17 H=2.D0*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.D0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.D0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
      END
      SUBROUTINE SIMPS2_Wdec(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C =====================================================================
C SIMPS
C A1,B1 -THE LIMITS OF INTEGRATION
C H1    -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN  AEPS1=REPS1=0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(7),P(5)
      external funct
      H=DSIGN(H1,B1-A1)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(X)/3.D0
    4 X0=X
      IF((X0+4.D0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.D0
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(X)/3.D0
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*ABS(F(K))
      DI1=(F(1)+4.D0*F(3)+F(5))*2.D0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.D0) 17,14,14
   17 H=2.D0*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.D0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.D0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
      END
      SUBROUTINE SIMPS3_Wdec(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C =====================================================================
C SIMPS
C A1,B1 -THE LIMITS OF INTEGRATION
C H1    -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN  AEPS1=REPS1=0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(7),P(5)
      external funct
      H=DSIGN(H1,B1-A1)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(X)/3.D0
    4 X0=X
      IF((X0+4.D0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.D0
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(X)/3.D0
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*ABS(F(K))
      DI1=(F(1)+4.D0*F(3)+F(5))*2.D0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.D0) 17,14,14
   17 H=2.D0*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.D0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.D0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
      END
c end of packege hard ========================================


