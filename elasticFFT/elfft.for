      program fft
c
      INCLUDE 'elfft.dim'
c
      integer p,q
c
      dimension data(2*npts1*npts2*npts3),nn(3),nn2(2)
      dimension xk(3)
      dimension a(3,3),g1(3,3,3,3)
c
      dimension delta(6,npts1,npts2,npts3)
      dimension deltaim(6,npts1,npts2,npts3)
      dimension aux6(6),aux33(3,3)
      dimension aux66(6,6),aux3333(3,3,3,3)
      dimension minv1(3),minv2(3)
c
      dimension d6(6),d6im(6),d33(3,3),d33im(3,3),fbar(3,3)
c
      dimension daux(6),dloc(3,3)
      dimension saux(6),sloc(3,3)
c
      dimension velmax(3)
      pi=4.*atan(1.)
c
      open(21,file='err.out',status='unknown')
cw      open(25,file='conv.out',status='unknown')
c
      nn(1)=npts1
      nn(2)=npts2
      nn(3)=npts3
c
      nn2(1)=npts1
      nn2(2)=npts2
c
      prodnn=float(nn(1)*nn(2)*nn(3))
      wgt=1./prodnn
c
      ur0=0
      open(ur0,file='srihari.in',status='old')
      UR1=1      ! FILECRYS
      UR2=2      ! FILETEXT
c
      call input
c
      write(21,'(a97)') ' IT  EERR       SERR       ETEN                      
     #                                              STEN'
      do 777 i=1,npts1
      do 777 j=1,npts2
      do 777 k=1,npts3
c
c     INITIALIZE D~,SG
c
      do ii=1,6
      dtilde(ii,i,j,k)=0.
      enddo
c
      do ii=1,6
c
cw      if(jphase(i,j,k).eq.0) then
      if(jphase(i,j,k).eq.2) then
c
         sg(ii,i,j,k)=0.
c
      else
c
         do jj=1,6
          sg(ii,i,j,k)=sg(ii,i,j,k)+
     #          cloc(ii,jj,i,j,k)*(dbar(jj)+dtilde(jj,i,j,k))
         enddo
c
      endif
c
         sbar(ii)=sbar(ii)+sg(ii,i,j,k)*wgt
c
      enddo

c      write(*,*) cloc(:,:,i,j,k)
c
777   continue
c
    
      call chg_basis(sbar,stens,aux66,aux3333,1)
      sref=stens(ictrl1,ictrl2)

      do 3000 imicro=1,nsteps

csshcw      write(*,*) 'STEP = ',imicro
cw      write(21,*) 'STEP = ',imicro
cw      write(25,*) 'STEP = ',imicro

      iter=0
c
      err2mod=2*error
      do while(iter.lt.itmax.and.err2mod.gt.error)
      iter=iter+1
c
cssh      write(*,*)'ITER = ',iter
c
cssh      write(*,*) 
cssh     #'DIRECT FFT OF POLARIZATION AND LAGRANGE MULTIPLIER FIELDS'
csshcw      write(*,*)
c
      do 300 ii=1,6
c
      k1=0
      do 5 k=1,npts3
c
cw      write(*,'(1H+,a,i2,2(a,i4))') 
cw     #    'DELTA - COMPONENT',ii,'  -   Z = ',k,'   OUT OF ',npts
c
      do 5 j=1,npts2
      do 5 i=1,npts1
      delta(ii,i,j,k)=sg(ii,i,j,k)
c
      do jj=1,6
      delta(ii,i,j,k)=delta(ii,i,j,k)-xlsec(ii,jj)*dtilde(jj,i,j,k)
      enddo
c
      k1=k1+1
      data(k1)=delta(ii,i,j,k)
      k1=k1+1
      data(k1)=0.
5     continue 
c
      if(npts3.gt.1) then
      call fourn(data,nn,3,1)
      else
      call fourn(data,nn2,2,1)
      endif
c
      k1=0
      do 6 kzz=1,npts3
      do 6 kyy=1,npts2
      do 6 kxx=1,npts1
c
      k1=k1+1
      delta(ii,kxx,kyy,kzz)=data(k1)
c
      k1=k1+1
      deltaim(ii,kxx,kyy,kzz)=data(k1)
c
6     continue
c
300   continue
c
cssh      write(*,*) 'CALCULATING G^pq,ij : TAU^ij ...'
csshcw      write(*,*)
c
      count=0;
      do 1 kzz=1,npts3
c
cw      write(*,'(1H+,a,i4,a,i4)') ' Z = ',kzz,'   OUT OF ',npts
c
      do 1 kyy=1,npts2
      do 1 kxx=1,npts1
c     
      count=count+1
      do 46 i=1,6
      d6(i)=delta(i,kxx,kyy,kzz)
      d6im(i)=deltaim(i,kxx,kyy,kzz)
46    continue
      if(iter .eq.1) then
      write(*,*) d6(:)
      write(*,*) d6im(:)
      endif

      call chg_basis(d6,d33,aux66,aux3333,1)
      call chg_basis(d6im,d33im,aux66,aux3333,1)
c

      if(kxx.le.npts1/2) kx=kxx-1
      if(kxx.gt.npts1/2) kx=kxx-npts1-1
c
      if(kyy.le.npts2/2) ky=kyy-1
      if(kyy.gt.npts2/2) ky=kyy-npts2-1
c
      if(kzz.le.npts3/2) kz=kzz-1
      if(kzz.gt.npts3/2) kz=kzz-npts3-1
c
      xk(1)=kx/(delt(1)*nn(1))
      xk(2)=ky/(delt(2)*nn(2))
c
      if(npts3.gt.1) then
      xk(3)=kz/(delt(3)*nn(3))
      else
      xk(3)=0.
      endif
c
      xknorm=sqrt(xk(1)**2+xk(2)**2+xk(3)**2)
c
      if (xknorm.ne.0.) then
      do i=1,3
c      xk2(i)=xk(i)/(xknorm*xknorm)
      xk(i)=xk(i)/xknorm
      enddo
      endif
c
      do 2 i=1,3
      do 2 k=1,3
      a(i,k)=0.
      do 2 j=1,3
      do 2 l=1,3
      a(i,k)=a(i,k)+xlsec33(i,j,k,l)*xk(j)*xk(l)
2     continue
c
      call minv(a,3,det,minv1,minv2)
c
c      if(det.eq.0) then
csshc      write(*,*) kx,ky,kz,'  --> SINGULAR SYSTEM'
c      stop
c      pause
c      endif
c
      do 3 p=1,3
      do 3 q=1,3
      do 3 i=1,3
      do 3 j=1,3
      g1(p,q,i,j)=-a(p,i)*xk(q)*xk(j)
3     continue
c
c      if (iter .eq. 1) then
csshc      write(*,*) count
c      do 1900 p=1,3
c      do 1900 q=1,3
csshc       write(*,*) g1(p,q,:,:)
c1900     continue
c      endif


c      if(iter .eq.1) then
c      write(*,*) kxx,kyy,kzz
c      write(*,*) d33
c      write(*,*) d33im
c      do 1900 p=1,3
c      do 1900 q=1,3
c      do i=1,3
c       write(*,*) (g1(p,q,i,j), j=1,3)
c      enddo
c      write(*,*)
c1900     continue
c
c      endif

      do 4 i=1,3
      do 4 j=1,3
c
      velgrad(i,j,kxx,kyy,kzz)=0.
      velgradim(i,j,kxx,kyy,kzz)=0.


c
      if(kx.eq.0.and.ky.eq.0.and.kz.eq.0.) goto 4
c
      do k=1,3
      do l=1,3
      velgrad(i,j,kxx,kyy,kzz)=
     #    velgrad(i,j,kxx,kyy,kzz)+g1(i,j,k,l)*d33(k,l)
      velgradim(i,j,kxx,kyy,kzz)=
     #    velgradim(i,j,kxx,kyy,kzz)+g1(i,j,k,l)*d33im(k,l)
      enddo
      enddo
c
4     continue 
c

c      if (iter .eq. 1) then
c     write(*,*) kxx,kyy,kzz
c      write(*,*) velgrad(:,:,kxx,kyy,kzz)
c      write(*,*) velgradim(:,:,kxx,kyy,kzz)
c      endif

1     continue
c
cc      call equilibrium(snormfft,snormfftim,nn,err2mod)
c
cssh      write(*,*) 'INVERSE FFT TO GET STRAIN FLUCTUATION FIELD'
c
      do 51 m=1,3
      do 51 n=1,3
c
      k1=0
c
      do 50 k=1,npts3
      do 50 j=1,npts2
      do 50 i=1,npts1
c
      k1=k1+1
      data(k1)=velgrad(m,n,i,j,k)
c
      k1=k1+1
      data(k1)=velgradim(m,n,i,j,k)
c
50     continue
c
      if(npts3.gt.1) then
      call fourn(data,nn,3,-1)
      else
      call fourn(data,nn2,2,-1)
      endif
c
      do i=1,2*npts1*npts2*npts3
      data(i)=data(i)/prodnn
      enddo
c
      k1=0
      do 16 kzz=1,npts3
      do 16 kyy=1,npts2
      do 16 kxx=1,npts1
      k1=k1+1
csshc      write(*,*) 'REAL PART =',kxx,kyy,kzz,data(k1)

      velgrad(m,n,kxx,kyy,kzz)=data(k1)
      k1=k1+1
csshc      write(*,*) 'IMAGINARY PART =',kxx,kyy,kzz,data(k1)
c      pause
16    continue
c
51    continue
c
c     DTILDE=SYM(VELGRAD)
c
      count=0
      do 167 kzz=1,npts3
      do 167 kyy=1,npts2
      do 167 kxx=1,npts1
c 
c      if (iter .eq. 1) then
c      write(*,*) kxx,kyy,kzz
c      write(*,*) velgrad(:,:,kxx,kyy,kzz)
c      endif

      count=count+1
      do ii=1,3
      do jj=1,3
      aux33(ii,jj)=(velgrad(ii,jj,kxx,kyy,kzz)+
     #              velgrad(jj,ii,kxx,kyy,kzz))/2.

      enddo
      enddo
c
      call chg_basis(aux6,aux33,aux66,aux3333,2)
c
c      write(*,*) aux6

      do m=1,6
      dtilde(m,kxx,kyy,kzz)=aux6(m)
      enddo

c      if(iter .eq. 1) then
csshc      write(*,*) count
csshc      write(*,*) dtilde(:,kxx,kyy,kzz)
c      endif
c
167   continue
c
cssh       write(*,*) 'UPDATE LAGRANGE MULTIPLIER FIELD'
       call augm_lagr
c
      do ii=1,6
      sbar(ii)=0.
      enddo
c
      do 757 i=1,npts1
      do 757 j=1,npts2
      do 757 k=1,npts3
c
c      if(iter .eq. 50) then
c      write(*,*) SG(:,I,J,K)
c      endif

      do ii=1,6
      sbar(ii)=sbar(ii)+sg(ii,i,j,k)*wgt
      enddo
c
757   continue
c
      call chg_basis(sbar,stens,aux66,aux3333,1)
      sref=stens(ictrl1,ictrl2)
c
cssh       write(*,*) 'STRESS FIELD ERROR =',errs/sref
cssh       write(*,*) 'STRAIN FIELD ERROR =',errd/dref
c
cw      write(21,101) iter,errd/dref,errs/sref,dref,sref
      write(21,101) iter,errd/dref,errs/sref,
     #       dsim(1,1),dsim(2,2),dsim(3,3),
     #       dsim(2,3),dsim(1,3),dsim(1,2),
     #       stens(1,1),stens(2,2),stens(3,3),
     #       stens(2,3),stens(1,3),stens(1,2)
101   format(i3,14(1x,e10.3))
c
c     ENDDO ... WHILE
c
      enddo
c      write(*,*) xlsec
cw
cw      IF(IUPDATE.EQ.1) THEN
cwc
cwc     VELMAX
cwc
cw      velmax(1)=dsim(1,1)*delt(1)*(npts1-1)
cw      velmax(2)=dsim(2,2)*delt(2)*(npts2-1)
cw      velmax(3)=dsim(3,3)*delt(3)*(npts3-1)
cwc
cwc     UPDATE DELT
cwc
cw      delt(1)=(delt(1)*(npts1-1)+velmax(1)*tdot)/(npts1-1)
cw      delt(2)=(delt(2)*(npts2-1)+velmax(2)*tdot)/(npts2-1)
cw      if(npts3.gt.1) then
cw      delt(3)=(delt(3)*(npts3-1)+velmax(3)*tdot)/(npts3-1)
cw      endif
cw
cw      ENDIF
cw
3000  CONTINUE
c
c
c     FIELDS.OUT
c
      open(24,file='fields.out',status='unknown')
c
      zero=0.
      one=1.
c
      do i=1,3
      do j=1,3
      fbar(i,j)=0.
      enddo
      enddo
c
      do i=1,3
      fbar(i,i)=delt(i)
      enddo
c
      write(24,'(i10)') npts1
      write(24,'(i10)') npts1*npts2*npts3
      write(24,111)((fbar(i,j),j=1,3),i=1,3)
      write(24,'(a93)') '   X   Y   Z  PH   GR  ELOC                      
     #                                          SLOC'
111   format(9f7.3) 
c
cw      ig=0
      do k=1,npts3
      do j=1,npts2
      do i=1,npts1
cw
cw      ig=ig+1
cw      do ii=1,3
cw      do jj=1,3
cw     aux33(ii,jj)=ag(jj,ii,i,j,k)
cw      enddo
cw      enddo
cwc
cw      call euler(1,ph,th,om,aux33)
cw
      do jj=1,6
      daux(jj)=dbar(jj)+dtilde(jj,i,j,k)
      saux(jj)=sg(jj,i,j,k)
      enddo
c
      call chg_basis(daux,dloc,aux66,aux3333,1)
      call chg_basis(saux,sloc,aux66,aux3333,1)
c
      write(24,172) i,j,k,jphase(i,j,k),jgrain(i,j,k),
     #       dloc(1,1),dloc(2,2),dloc(3,3),
     #       dloc(2,3),dloc(1,3),dloc(1,2),
     #       sloc(1,1),sloc(2,2),sloc(3,3),
     #       sloc(2,3),sloc(1,3),sloc(1,2)
c
172   format(4i4,i5,12e11.3)
c
      enddo
      enddo
      enddo
c
      end
c     
      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      REAL data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
     *k2,n,nprev,nrem,ntot
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then

            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif

          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr

                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END
c
      SUBROUTINE VOIGT(C2,C4,IOPT)
C
C *** TRANSFORMS SECOND ORDER MATRIX C2 INTO FOURTH ORDER TENSOR C4 IF
C *** IOPT=1 AND VICEVERSA IF IOPT=2. IF IOPT=3,TRANSFORMS WITH INV.FACT.
C *** IOPT=4 FOR GO FROM 6x6 TO 3x3x3x3 WITH Aijkl ANTISYMMETRY
C
      DIMENSION C2(6,6),C4(3,3,3,3),IJV(6,2),F(6,6)
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
C
      IF(IOPT.EQ.1) THEN
      DO 10 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 10 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
      C4(I1,I2,J1,J2)=C2(I,J)
      C4(I2,I1,J1,J2)=C2(I,J)
      C4(I1,I2,J2,J1)=C2(I,J)
   10 C4(I2,I1,J2,J1)=C2(I,J)
      ENDIF
C
      IF(IOPT.EQ.2) THEN
      DO 20 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 20 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
   20 C2(I,J)=C4(I1,I2,J1,J2)
      ENDIF
c      
      IF(IOPT.EQ.3) THEN
      DO 9 I=1,6
      DO 9 J=1,6
      F(I,J)=1.
      IF(I.GT.3) F(I,J)=0.5
      IF(J.GT.3) F(I,J)=0.5*F(I,J)
9     CONTINUE
C
      DO 101 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 101 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
      C4(I1,I2,J1,J2)=F(I,J)*C2(I,J)
      C4(I2,I1,J1,J2)=F(I,J)*C2(I,J)
      C4(I1,I2,J2,J1)=F(I,J)*C2(I,J)
101   C4(I2,I1,J2,J1)=F(I,J)*C2(I,J)
      ENDIF
C
      IF(IOPT.EQ.4) THEN
      DO 17 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 17 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
      IF(I.LE.3) THEN
      C4(I1,I2,J1,J2)=C2(I,J)
      C4(I2,I1,J1,J2)=C2(I,J)
      C4(I1,I2,J2,J1)=C2(I,J)
      C4(I2,I1,J2,J1)=C2(I,J)
      ELSE
      C4(I1,I2,J1,J2)=C2(I,J)
      C4(I2,I1,J1,J2)=-C2(I,J)
      C4(I1,I2,J2,J1)=C2(I,J)
      C4(I2,I1,J2,J1)=-C2(I,J)
      ENDIF
17    CONTINUE
      ENDIF
      RETURN
      END
c
      subroutine euler(iopt,ph,th,tm,a)
      dimension a(3,3)
      pi=4.*atan(1.d0)
c
c     CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
c     MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
c     A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.
c     ph,th,om ARE THE EULER ANGLES OF ca REFERRED TO sa.
c
      if(iopt.eq.1) then
        th=acos(a(3,3))
        if(abs(a(3,3)).ge.0.9999) then
          tm=0.
          ph=atan2(a(1,2),a(1,1))
        else
          sth=sin(th)
          tm=atan2(a(1,3)/sth,a(2,3)/sth)
          ph=atan2(a(3,1)/sth,-a(3,2)/sth)
        endif
        th=th*180./pi
        ph=ph*180./pi
        tm=tm*180./pi
      else if(iopt.eq.2) then
        sph=sin(ph)
        cph=cos(ph)
        sth=sin(th)
        cth=cos(th)
        stm=sin(tm)
        ctm=cos(tm)
        a(1,1)=ctm*cph-sph*stm*cth
        a(2,1)=-stm*cph-sph*ctm*cth
        a(3,1)=sph*sth
        a(1,2)=ctm*sph+cph*stm*cth
        a(2,2)=-sph*stm+cph*ctm*cth
        a(3,2)=-sth*cph
        a(1,3)=sth*stm
        a(2,3)=ctm*sth
        a(3,3)=cth
      endif

      return
      end
c
C ********************************************************************
C     SUBROUTINE INPUT      --->      VERSION 31/jan/99
C
C     READS CHARACTERISTICS OF THE RUN: # OF PHASES, NAMES OF INPUT FILES,
C     DEFORMATION TO BE IMPOSED, CONVERGENCE PARAMETERS, ETC.
C     READS SINGLE CRYSTAL PROPERTIES: DEFORMATION MODES, CRSS, HARDENING
C     READS CRYSTAL AND MORPHOLOGIC TEXTURES.
C     INITIALIZES ARRAYS REQUIRED TO RUN VPSC.
C     OPENS AND CLOSES INPUT FILES.   OPENS OUTPUT FILES.
C
C     MODIFIED 21/07/98 by CNT:
C     INITIALIZATION RELATED TO 'ELEMENTS' IS DONE INSIDE A SINGLE BLOCK.
C *****************************************************************************
C
      SUBROUTINE INPUT

      INCLUDE 'elfft.dim'

      dimension aux66(6,6),aux3333(3,3,3,3)
      DIMENSION IJV(6,2)
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/

C *********   INITIALIZATION BLOCK   ***************************
C
      PI=4.*ATAN(1.)

C     CALCULATES TENSORS OF THE SYMMETRIC BASIS 'B(3,3,6)'

      CALL CHG_BASIS(DUM1,DUM2,DUM3,DUM4,0)
C
      READ(UR0,'(a)') prosa
      READ(UR0,'(a)') filetext
      READ(UR0,'(a)') prosa
      READ(UR0,'(a)') filecrys
c
      OPEN (unit=UR1,file=filecrys,status='old')
        call data_crystal
      CLOSE(unit=UR1)
c
      OPEN(unit=UR2,file=filetext,status='old')
        call data_grain
      CLOSE(unit=UR2)
C
      READ(UR0,'(A)') PROSA
c
c     RVE DIMENSIONS
c
      READ(UR0,*) DELT

C ****************************************************************************
C     READ BOUNDARY CONDITIONS ON OVERALL STRAIN
C ****************************************************************************

      READ(UR0,'(A)') PROSA
C
      iudot=1
c
      DO I=1,3
        READ(UR0,*) (UDOT(I,J),J=1,3)
      ENDDO
C
c     SYMMETRIC STRAIN, ANTISYMMETRIC ROTATION TENSORS
c     AND INDICES OF IMPOSED COMPONENTS
c
      do i=1,3
      do j=1,3
        dsim(i,j)=(udot(i,j)+udot(j,i))/2.
        tomtot(i,j)=(udot(i,j)-udot(j,i))/2.
      enddo
      enddo
c
      do i=1,3
        idsim(i)=iudot(i,i)
      enddo
c
      idsim(4)=0
      if(iudot(2,3).eq.1.and.iudot(3,2).eq.1) idsim(4)=1
      idsim(5)=0
      if(iudot(1,3).eq.1.and.iudot(3,1).eq.1) idsim(5)=1
      idsim(6)=0
      if(iudot(1,2).eq.1.and.iudot(2,1).eq.1) idsim(6)=1
c
c     WRITES STRAIN DSIM(I,J) IN b-BASIS AS A 5-DIM VECTOR DBAR(K)
c
      call chg_basis(dbar,dsim,aux66,aux3333,2)
c
      READ(UR0,'(A)') PROSA
cw      READ(UR0,*) DUMMY
c
      READ(UR0,*) ICTRL
c
      if(ictrl.le.3) then
      ictrl1=ictrl
      ictrl2=ictrl
      else if(ictrl.eq.4) then
      ictrl1=2
      ictrl2=3
      else if(ictrl.eq.5) then
      ictrl1=1
      ictrl2=3
      else if(ictrl.eq.6) then
      ictrl1=1
      ictrl2=2
      endif
c
      dref=dsim(ictrl1,ictrl2)
c
cw      if(ictrl.eq.0)  evmincr=dummy
cw      if(ictrl.gt.0)  eijincr=dummy
c
      READ(UR0,*) NSTEPS
cw      nsteps=1

      READ(UR0,*) ERROR
      READ(UR0,*) ITMAX
c
      RETURN
      END
c
C *****************************************************************************
C     SUBROUTINE DATA_CRYSTAL        --->      VERSION 03/FEB/2000
C *****************************************************************************

      SUBROUTINE DATA_CRYSTAL
c
      INCLUDE 'elfft.dim'
c
	dimension dde(3,3),xid4(3,3,3,3)
cw	dimension aux6(6),aux33(3,3)
c
c       UNITARY TENSORS
c
	do i=1,3
	do j=1,3
	dde(i,j)=0.d0
	if(i.eq.j) dde(i,j)=1.d0
	enddo
	enddo
c
	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
	xid4(i,j,k,l)=(dde(i,k)*dde(j,l)+dde(i,l)*dde(j,k))/2.d0
	enddo
	enddo
	enddo
	enddo
c
cw	read(UR1,*) fact
        fact=0.5
	read(UR1,*) iso
c
	if(iso.eq.0) then
c
	do i=1,6
	read(UR1,*)(cmat(i,j),j=1,6)
	enddo
	call voigt(cmat,cmat33,1)
c
	else
c
cw	read(UR1,*) tmu,tnu
	read(UR1,*) young,tnu
        tmu=young/(2.*(1.+tnu))
	tla=2.d0*tmu*tnu/(1.d0-2.d0*tnu)
c       bulk=2.d0*tmu*(1+tnu)/3.d0/(1.d0-2d0*tnu)
c
c	rho=3.d0*(1.d0-2.d0*tnu)/(1.d0+tnu)
c
	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
	cmat33(i,j,k,l)=tla*dde(i,j)*dde(k,l)+2.d0*tmu*xid4(i,j,k,l)
	enddo
	enddo
	enddo
	enddo
c
     	endif
c
cw      call chg_basis(aux6,aux33,cmat,cmat33,4)
c
      RETURN
      END
C
      SUBROUTINE DATA_GRAIN

      INCLUDE 'elfft.dim'
c
      DIMENSION AA(3,3)
      DIMENSION CAUX33(3,3,3,3),CAUX(6,6)
      dimension saux(6,6),taux(6,6)
      dimension aux6(6),aux33(3,3)
      dimension minv1(6),minv2(6)
c
      do i=1,6
      do j=1,6
      xlsec(i,j)=0.
      enddo
      enddo
c
      do kkk=1,npts1*npts2*npts3
c
      READ(UR2,*) ph,th,om,ii,jj,kk,jgr,jph
c
      jgrain(ii,jj,kk)=jgr
      jphase(ii,jj,kk)=jph
C
C     CALCULATES THE TRANSFORMATION 
C     MATRIX AA, WHICH TRANSFORMS FROM SAMPLE TO CRYSTAL.
C     AG TRANSFORMS FROM CRYSTAL TO SAMPLE.

        CALL EULER(2,ph*pi/180,th*pi/180,om*pi/180,aa)
c
cw        DO J=1,3
cw        DO K=1,3
cw          AG(J,K,ii,jj,kk)=AA(K,J)
cw        ENDDO
cw        ENDDO
c
      do 1 i1=1,3
      do 1 j1=1,3
      do 1 k1=1,3
      do 1 l1=1,3
      dum=0.
      do 2 i2=1,3
      do 2 j2=1,3
      do 2 k2=1,3
      do 2 l2=1,3
      dum=dum+
     #    aa(i2,i1)*aa(j2,j1)*aa(k2,k1)*aa(l2,l1)*cmat33(i2,j2,k2,l2)
2     continue
      caux33(i1,j1,k1,l1)=dum
1     continue
c
      call chg_basis(aux6,aux33,caux,caux33,4)
c
      do i=1,6
      do j=1,6
      cloc(i,j,ii,jj,kk)=caux(i,j)
      xlsec(i,j)=xlsec(i,j)+caux(i,j)*wgt
      enddo
      enddo
cwritecaux      write(*,*) caux
c
      enddo
c
      call chg_basis(aux6,aux33,xlsec,xlsec33,3)      
C
      do ii=1,npts1
      do jj=1,npts2
      do kk=1,npts3
c
      do i=1,6
      do j=1,6
      saux(i,j)=cloc(i,j,ii,jj,kk)
      enddo
      enddo
c
      call minv(saux,6,det,minv1,minv2)
c 
      do i=1,6
      do j=1,6
      dum=0.
      do k=1,6
      dum=dum+xlsec(i,k)*saux(k,j)
      enddo
      taux(i,j)=(i/j)*(j/i)+dum
      enddo
      enddo
c
      call minv(taux,6,det,minv1,minv2)
c
      do i=1,6
      do j=1,6
      dum=0.
      do k=1,6
      dum=dum+saux(i,k)*taux(k,j)
      enddo
      fsloc(i,j,ii,jj,kk)=dum
      enddo
      enddo

csshc      write(*,*) fsloc(:,:,ii,jj,kk)
c
      enddo
      enddo
      enddo
c
csshc      write(*,*) xlsec
      RETURN
      END

C ************************************************************************
C
C     SUBROUTINE CHG_BASIS_OLD    --->   VERSION 06/02/98
C
C     (modif. 10/JAN/97 - KDIM version - R.L.)
C     (modif. 06/FEB/98 - same convention as SELFPOLY - C.N.T.)
C     (modif. 16/JUN/99 - same convention as Maudlin  - C.N.T.)
C
C     DEFINES A BASIS OF 6 SECOND ORDER TENSORS B(N) IF IOPT=0.
C     CALCULATES THE SECOND ORDER TENSOR 'C2' AS AN EXPANSION IN
C     TERMS OF ITS COMPONENTS CE2(N) AND THE BASIS TENSORS B(N) IF
C     IOPT=1.
C     CALCULATES THE COMPONENTS OF C2 AS A VECTOR CE2(6) IF IOPT=2.
C     CALCULATES THE FOURTH ORDER TENSOR 'C4' AS AN EXPANSION IN TERMS
C     OF ITS COMPONENTS CE4(N,M) AND THE BASIS TENSORS B(N) IF IOPT=3.
C     CALCULATES THE COMPONENTS OF C4 AS A MATRIX CE4(6,6) IF IOPT=4.
C **************************************************************************
C
      SUBROUTINE CHG_BASIS(CE2,C2,CE4,C4,IOPT)
C
      PARAMETER (KDIM=6)
c      PARAMETER (SQR2=1.41421356237309   )
      PARAMETER (RSQ2=0.70710678118654744)
      PARAMETER (RSQ3=0.57735026918962584)
      PARAMETER (RSQ6=0.40824829046386304)
C
      DIMENSION CE2(KDIM),C2(3,3),CE4(KDIM,KDIM),C4(3,3,3,3)
C
C     DIMENSION B(3,3,6)
C     DATA B /RSQ6,0,   0,   0,   RSQ6,0,   0,   0,  -2*RSQ6,
C    #        RSQ2,0,   0,   0,  -RSQ2,0,   0,   0,   0,
C    #        0,   0,   0,   0,   0,   RSQ2,0,   RSQ2,0,
C    #        0,   0,   RSQ2,0,   0,   0,   RSQ2,0,   0,
C    #        0,   RSQ2,0,   RSQ2,0,   0,   0,   0,   0,
C    #        RSQ3,0,   0,   0,   RSQ3,0,   0,   0,   RSQ3/

      COMMON/BASIS/ B(3,3,6)

      IF(IOPT.EQ.0) THEN
C *** CALCULATES BASIS TENSORS B(N)

        DO I=1,3
          DO J=1,3
            DO N=1,6
              B(I,J,N)=0.0
            ENDDO
          ENDDO
        ENDDO

        B(1,1,2)=-RSQ6
        B(2,2,2)=-RSQ6
        B(3,3,2)= 2.D0*RSQ6

        B(1,1,1)=-RSQ2
        B(2,2,1)= RSQ2

        B(2,3,3)=RSQ2
        B(3,2,3)=RSQ2

        B(1,3,4)=RSQ2
        B(3,1,4)=RSQ2

        B(1,2,5)=RSQ2
        B(2,1,5)=RSQ2

        B(1,1,6)=RSQ3
        B(2,2,6)=RSQ3
        B(3,3,6)=RSQ3

      ENDIF
C
      IF(IOPT.EQ.1) THEN
C *** CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
      DO 40 I=1,3
      DO 40 J=1,3
      C2(I,J)=0.0
      DO 40 N=1,KDIM
   40 C2(I,J)=C2(I,J)+CE2(N)*B(I,J,N)
C
      ENDIF
C
      IF(IOPT.EQ.2) THEN
C *** CALCULATES KDIMx1 b-COMPONENTS VECTOR FROM SECOND ORDER TENSOR.
      DO 50 N=1,KDIM
      CE2(N)=0.0
      DO 50 I=1,3
      DO 50 J=1,3
   50 CE2(N)=CE2(N)+C2(I,J)*B(I,J,N)
C
      ENDIF
C
      IF(IOPT.EQ.3) THEN
C *** CALCULATES FOURTH ORDER TENSOR FROM b-COMPONENTS MATRIX.
      DO 20 I=1,3
      DO 20 J=1,3
      DO 20 K=1,3
      DO 20 L=1,3
      C4(I,J,K,L)=0.0
      DO 20 N=1,KDIM
      DO 20 M=1,KDIM
   20 C4(I,J,K,L)=C4(I,J,K,L)+CE4(N,M)*B(I,J,N)*B(K,L,M)
C
      ENDIF
C
      IF(IOPT.EQ.4) THEN
C *** CALCULATES KDIMxKDIM b-COMPONENTS MATRIX FROM FOURTH ORDER TENSOR.
      DO 30 N=1,KDIM
      DO 30 M=1,KDIM
      CE4(N,M)=0.0
      DO 30 I=1,3
      DO 30 J=1,3
      DO 30 K=1,3
      DO 30 L=1,3
   30 CE4(N,M)=CE4(N,M)+C4(I,J,K,L)*B(I,J,N)*B(K,L,M)
C
      ENDIF
C
      RETURN
      END
c
      SUBROUTINE AUGM_LAGR

      INCLUDE 'elfft.dim'
c
      dimension x(6),dg(6),edot(6),dsg(6),ddg(6)
c
       ERRD=0.
       ERRS=0.
c
       do kk=1,npts3
c
cw       write(*,'(1H+,a,i4,a,i4)') ' Z = ',kk,'   OUT OF ',npts
c
       do jj=1,npts2
       do ii=1,npts1
c
       if(jphase(ii,jj,kk).eq.1) then
c
        DO I=1,6
        DG(I)=DBAR(I)+DTILDE(I,II,JJ,KK)
        ENDDO

c
        DO I=1,6
          X(I)=SG(I,II,JJ,KK)
          DO J=1,6
            X(I)=X(I)+XLSEC(I,J)*DG(J)
          ENDDO
        ENDDO

c
        DO I=1,6
          EDOT(I)=0
          DO J=1,6
            EDOT(I)=EDOT(I)+FSLOC(I,J,II,JJ,KK)*X(J)
          ENDDO
        ENDDO

c     UPDATE SG (LAGRANGE MULTIPLIER)

      DO I=1,6
       DDG(I)=DG(I)-EDOT(I)
       DSG(I)=0.
       DO J=1,6
        DSG(I)=DSG(I)+XLSEC(I,J)*(DG(J)-EDOT(J))
       ENDDO
       SG(I,II,JJ,KK)=SG(I,II,JJ,KK)+DSG(I)
      ENDDO
c
      ERRD=ERRD+TNORM(DDG,6,1)*WGT
      ERRS=ERRS+TNORM(DSG,6,1)*WGT
c
      endif      ! jphase endif
c
      ENDDO
      ENDDO
      ENDDO
c
      RETURN
      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     FUNCTION TMISMATCH   ---->   VERSION OF 27/DEC/98
C
C     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO NROWSxNCOLS MATRICES
C     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE
C     OF BOTH DATA.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION tmismatch (v1,v2,nrows,ncols)
C
COLD      DIMENSION v1(nrows,ncols),v2(nrows,ncols)
COLD      DIMENSION v_dif(6,6),v_ave(6,6)
      DIMENSION v1(36),v2(36)
      DIMENSION v_dif(36),v_ave(36)
C
      do i=1,nrows*ncols
        v_dif(i)=v1(i)-v2(i)
        v_ave(i)=0.5d0*(v1(i)+v2(i))
      enddo
      tmismatch=tnorm(v_dif,nrows,ncols)/tnorm(v_ave,nrows,ncols)
C
      RETURN
      END

C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     FUNCTION TNORM   ---->   VERSION OF 27/DEC/98
C
C     CALCULATES THE NORM OF A NROWSxNCOLS-MATRIX (NROWS,NCOLS =< 6)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION tnorm(v,nrows,ncols)
C
COLD      DIMENSION v(nrows,ncols)
      DIMENSION v(36)
C
      tnorm=0.d0
      do i=1,nrows*ncols
        tnorm=tnorm+v(i)*v(i)
      enddo
      tnorm=sqrt(tnorm)
C
      RETURN
      END
c
      SUBROUTINE MINV (A,N,D,L,M)
c
      DIMENSION A(*),L(*),M(*)
C
C       SEARCH FOR LARGEST ELEMENT
C
      D=1.d0
      NK=-N
      DO 180 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
      IF (ABS(BIGA)-ABS(A(IJ))) 10,20,20
   10 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
C
C       INTERCHANGE ROWS
C
      J=L(K)
      IF (J-K) 50,50,30
   30 KI=K-N
      DO 40 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   40 A(JI)=HOLD
C
C       INTERCHANGE COLUMNS
C
   50 I=M(K)
      IF (I-K) 80,80,60
   60 JP=N*(I-1)
      DO 70 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   70 A(JI)=HOLD
C
C       DIVIDE COLUMN BY MINUS PIVOT (BIGA)
C
   80 IF (ABS(BIGA).LT.1.e-10) THEN
   90 D=0.0
      RETURN
      ENDIF
  100 DO 120 I=1,N
      IF (I-K) 110,120,110
  110 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
  120 CONTINUE
C
C       REDUCE MATRIX
C
      DO 150 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 150 J=1,N
      IJ=IJ+N
      IF (I-K) 130,150,130
  130 IF (J-K) 140,150,140
  140 KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
  150 CONTINUE
C
C       DIVIDE ROW BY PIVOT
C
      KJ=K-N
      DO 170 J=1,N
      KJ=KJ+N
      IF (J-K) 160,170,160
  160 A(KJ)=A(KJ)/BIGA
  170 CONTINUE
C
C       PRODUCT OF PIVOTS
C
      D=D*BIGA
C
C       REPLACE PIVOT BY RECIPROCAL
C
      A(KK)=1.d0/BIGA
  180 CONTINUE
C
C       FINAL ROW AND COLUMN INTERCHANGE
C
      K=N
  190 K=(K-1)
      IF (K) 260,260,200
  200 I=L(K)
      IF (I-K) 230,230,210
  210 JQ=N*(K-1)
      JR=N*(I-1)
      DO 220 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  220 A(JI)=HOLD
  230 J=M(K)
      IF (J-K) 190,190,240
  240 KI=K-N
      DO 250 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  250 A(JI)=HOLD
      GO TO 190
  260 RETURN
      END
