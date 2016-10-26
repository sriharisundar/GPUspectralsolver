! VERSION 7  with FFTW lib
! Anand Kanjarla
! (i) voce type hardening -- mixed implicit explicit procedure
! (ii) texture evoltuion is corrected - compared with experimental texture  
!
! adding differenct criteria for convergene of stress
! (i) divergence in fourier space
!  
! 
!**************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! this module defines in a proper way the variables associated to the different precision kinds    
   module precision_kinds    
     implicit none
     integer ( kind = kind (1) ) , parameter :: i2b = kind ( 1 ) !> simple precision integer
     integer ( kind = i2b ) , parameter :: dp = kind ( 0.0d0 ) !> double precision real
     integer ( kind = i2b ) , parameter :: sp = kind ( 0.0 ) !> simple precision real    
     integer ( kind = i2b ) , parameter :: i4b = 2_i2b * i2b !> double precision integer
   end module precision_kinds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***************************************************************
 module fftw3_f90
  integer ( kind = 4 ), parameter :: fftw_r2hc = 0
  integer ( kind = 4 ), parameter :: fftw_hc2r = 1
  integer ( kind = 4 ), parameter :: fftw_dht = 2
  integer ( kind = 4 ), parameter :: fftw_redft00 = 3
  integer ( kind = 4 ), parameter :: fftw_redft01 = 4
  integer ( kind = 4 ), parameter :: fftw_redft10 = 5
  integer ( kind = 4 ), parameter :: fftw_redft11 = 6
  integer ( kind = 4 ), parameter :: fftw_rodft00 = 7
  integer ( kind = 4 ), parameter :: fftw_rodft01 = 8
  integer ( kind = 4 ), parameter :: fftw_rodft10 = 9
  integer ( kind = 4 ), parameter :: fftw_rodft11 = 10
  integer ( kind = 4 ), parameter :: fftw_forward = -1
  integer ( kind = 4 ), parameter :: fftw_backward = +1
  integer ( kind = 4 ), parameter :: fftw_measure = 0
  integer ( kind = 4 ), parameter :: fftw_destroy_input = 1
  integer ( kind = 4 ), parameter :: fftw_unaligned = 2
  integer ( kind = 4 ), parameter :: fftw_conserve_memory = 4
  integer ( kind = 4 ), parameter :: fftw_exhaustive = 8
  integer ( kind = 4 ), parameter :: fftw_preserve_input = 16
  integer ( kind = 4 ), parameter :: fftw_patient = 32
  integer ( kind = 4 ), parameter :: fftw_estimate = 64
  integer ( kind = 4 ), parameter :: fftw_estimate_patient = 128
  integer ( kind = 4 ), parameter :: fftw_believe_pcost = 256
  integer ( kind = 4 ), parameter :: fftw_dft_r2hc_icky = 512
  integer ( kind = 4 ), parameter :: fftw_nonthreaded_icky = 1024
  integer ( kind = 4 ), parameter :: fftw_no_buffering = 2048
  integer ( kind = 4 ), parameter :: fftw_no_indirect_op = 4096
  integer ( kind = 4 ), parameter :: fftw_allow_large_generic = 8192
  integer ( kind = 4 ), parameter :: fftw_no_rank_splits = 16384
  integer ( kind = 4 ), parameter :: fftw_no_vrank_splits = 32768
  integer ( kind = 4 ), parameter :: fftw_no_vrecurse = 65536
  integer ( kind = 4 ), parameter :: fftw_no_simd = 131072

end module 
!*******************************************************************
MODULE GLOBAL
use precision_kinds , only : i2b , i4b , dp
!INTEGER(kind =i2b) :: NPTS1,NPTS2,NPTS3
integer :: npts1,npts2,npts3
INTEGER, PARAMETER :: NPHMX=1,NMODMX=3,NTWMMX=1,NSYSMX=18
!      PARAMETER(NPTS1=32,NPTS2=32,NPTS3=32)
!      PARAMETER(NPHMX=2)     ! MAXIMUM # OF PHASES
!      PARAMETER(NMODMX=3)     ! MAXIMUM # OF ACTIVE SL+TW MODES IN ANY PHASE
!      PARAMETER(NTWMMX=1)     ! MAXIMUM # OF ACTIVE TWIN MODES IN ANY PHASE
!      PARAMETER(NSYSMX=12)    ! MAXIMUM # OF ACTIVE SL+TW SYSTEMS IN ANY PHASE
CHARACTER*80 FILETEXT,FILECRYSPL,FILECRYSEL,FILEHIST,PROSA
CHARACTER*3  ICRYST(NPHMX)
!INTEGER   ::   UR0,UR1,UR2,UR3,UR4,UR5,UR6,UW1,UW2

! TESTCOND      
DOUBLE PRECISION :: UDOT(3,3),DSIM(3,3),SCAUCHY(3,3),SDEVIAT(3,3),TOMTOT(3,3),DBAR5(5)
DOUBLE PRECISION :: DBAR6(6),DELT(3),DELTVOL3,TDOT, UDOT_IN(3,3)
DOUBLE PRECISION :: DEFGRADMACRO(3,3),DDEFGRADMACRO(3,3),scauav(3,3), DDEFGRADMACROACUM(3,3)
INTEGER :: IUDOT(3,3),IDSIM(6),ISCAU(6),ICTRL,NSTEPS
! DATACRYST
DOUBLE PRECISION :: DNCA(3,NSYSMX,NPHMX),DBCA(3,NSYSMX,NPHMX),SCHCA(5,NSYSMX,NPHMX)
DOUBLE PRECISION :: TAU(NSYSMX,3,NPHMX),HARD(NSYSMX,NSYSMX,NPHMX)
DOUBLE PRECISION :: THET(NSYSMX,1:2,NPHMX),GAMD0(NSYSMX,NPHMX)
INTEGER :: NRS(NSYSMX,NPHMX)
!DATAMODES
DOUBLE PRECISION :: TWSH(NSYSMX,NPHMX)
INTEGER :: NSM(NMODMX,NPHMX), NMODES(NPHMX),NSYST(NPHMX),NTWMOD(NPHMX),NTWSYS(NPHMX),ISECTW(NSYSMX,NPHMX)!,ICRYST(NPHMX)
!DATAPHASE
INTEGER :: NPH,NELEM,IPHBOT,IPHTOP,NGR
!FGRID
DOUBLE PRECISION :: WGT !,CRSS(NSYSMX,2,NPTS1,NPTS2,NPTS3),SG(3,3,NPTS1,NPTS2,NPTS3),EDOTP(3,3,NPTS1,NPTS2,NPTS3),EPT(3,3,NPTS1,NPTS2,NPTS3)
!REAL :: SCH(5,NSYSMX,NPTS1,NPTS2,NPTS3),AG(3,3,NPTS1,NPTS2,NPTS3),DEFGRAD(3,3,NPTS1,NPTS2,NPTS3),JPHASE(NPTS1,NPTS2,NPTS3),JGRAIN(NPTS1,NPTS2,NPTS3)
DOUBLE PRECISION , ALLOCATABLE :: CRSS(:,:,:,:,:),EDOTP(:,:,:,:,:),EPT(:,:,:,:,:),VELGRAD(:,:,:,:,:)
DOUBLE PRECISION , ALLOCATABLE :: AG(:,:,:,:,:),DEFGRAD(:,:,:,:,:)!,JPHASE(:,:,:),JGRAIN(:,:,:)
DOUBLE PRECISION , ALLOCATABLE :: RSST(:,:,:,:),SCH(:,:,:,:,:),DISGRAD_PRE(:,:,:,:,:),DISGRAD_CUR(:,:,:,:,:)
INTEGER, ALLOCATABLE :: JPHASE(:,:,:),JGRAIN(:,:,:)
REAL ( kind = dp ) , ALLOCATABLE, DIMENSION (:, :, : , : , : ) :: SG 
!DOUBLE PRECISION , ALLOCATABLE :: SG(:,:,:,:,:)
! IOFILES
!REAL :: FILETEXT,FILECRYSPL,FILECRYSEL,FILEHIST,PROSA
! IOUNITS
INTEGER :: UR0,UR1,UR2,UR3,UR4,UR5,UR6,UW1,UW2
! MISCEL
!REAL :: PI
INTEGER :: JRAN
! RUNCOND 
DOUBLE PRECISION :: ERROR,FACT2
INTEGER :: ITMAX, IRECOVER, ISAVE,IWTEX,IWDEQ,IFREQ
! PPC
DOUBLE PRECISION :: SVM,DVM,ERRS,ERRE,ERRE2, DIV_ERR
INTEGER :: IUPDATE,IUPHARD,IGAS(NPHMX),IWFIELDS,ILATTICE
! ELAS
DOUBLE PRECISION :: cc(3,3,3,3),c0(3,3,3,3),s0(3,3,3,3),c066(6,6) !cg66(6,6,NPTS1,NPTS2,NPTS3)
DOUBLE PRECISION, ALLOCATABLE :: cg66(:,:,:,:,:),GAMDOT(:,:,:,:), div_array(:,:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: div_array_im(:,:,:,:),div_REAL(:,:,:,:),PLAS_SLIP(:,:,:,:,:)
DOUBLE PRECISION :: tot_elas_ener, twin_frame(5)
!
! MISC ADDED BY ANAND 
!
DOUBLE PRECISION :: DIV_DOMAIN(3)
DOUBLE PRECISION,PARAMETER ::   Pi=3.141592653589793
DOUBLE PRECISION ::  PRODNN, EVM, evmp, dvmp ! PI,
INTEGER :: IMICRO
INTEGER :: ii,jj,kk,ll,iter,l !, K1
DOUBLE PRECISION :: COVERA
INTEGER :: NMODESX

DOUBLE PRECISION :: velmax(3), fbar(3,3)
DOUBLE PRECISION :: ph,th,om

DOUBLE PRECISION :: VOCE, FACT, EXPINI, EXPDEL, TAU0, TAU1,THET0,THET1,TINY
DOUBLE PRECISION , ALLOCATABLE  ::  GACUMGR(:,:,:),GACUMGRaux(:,:,:)
DOUBLE PRECISION , ALLOCATABLE :: TRIALTAU(:,:,:,:,:)


INTEGER :: modex,nsmx,nrsx(NSYSMX),isectwx, MODE(NSYSMX)
DOUBLE PRECISION :: gamd0x,twshx,tau0xf,tau0xb,tau1x,thet0x,thet1x
DOUBLE PRECISION :: HSELFX(12),HLATEX(12,12)

DOUBLE PRECISION :: JXRSINI,JXRSFIN,JXRSTEP
INTEGER :: JXRS,jph

DOUBLE PRECISION :: sinc, filter

DOUBLE PRECISION :: temp_shear(3,3)
!real :: shearstrain(3,3),temp_shear(3,3)
!real :: char_shear

! average values 

DOUBLE PRECISION :: stress_twin(6),stress_parent(6)
DOUBLE PRECISION :: elas_twin(6),elas_parent(6)
DOUBLE PRECISION:: plas_twin(6),plas_parent(6)
DOUBLE PRECISION :: el_ener_twin,el_ener_parent
DOUBLE PRECISION :: shear_twin,shear_parent
INTEGER :: vol_twin , vol_parent

INTEGER :: n_inc

CHARACTER(len=25) :: output_prefix
CHARACTER(len=25) :: CI, filename


END MODULE

!**************************************************************
!**************************************************************
!**************************************************************
PROGRAM FFTW_EIGEN_TWIN2a
use fftw3_f90
use GLOBAL
use precision_kinds , only : i2b , i4b , dp 

IMPLICIT NONE

!include '/usr/local/fftw-3.2.2/include/fftw3.f'
!include 'fftw3.f90'

integer p,q

INTEGER :: nn(3),nn2(2)
!DOUBLE PRECISION , ALLOCATABLE :: data(:)
DOUBLE PRECISION :: xk(3),xk2(3)
DOUBLE PRECISION :: a(3,3),g1(3,3,3,3)

DOUBLE PRECISION :: aux6(6),aux33(3,3)

DOUBLE PRECISION :: minv1(3),minv2(3)
!REAL (kind =dp) :: ddefgrad(3,3),ddefgradim(3,3)
DOUBLE PRECISION :: ddefgrad(3,3),ddefgradim(3,3)
DOUBLE PRECISION :: ddefgradav(3,3)

!REAL :: work(3,3,NPTS1,NPTS2,NPTS3)
!REAL :: workim(3,3,NPTS1,NPTS2,NPTS3)
!REAL ( KIND = dp ) , ALLOCATABLE, DIMENSION (:, :,: , : , : ) :: work,workim
DOUBLE PRECISION , ALLOCATABLE, DIMENSION(:,:,:,:,:) :: work, workim 
!DOUBLE PRECISION , ALLOCATABLE :: WORK(:,:,:,:,:),WORKIM(:,:,:,:,:)
DOUBLE PRECISION , ALLOCATABLE :: strain(:,:,:,:),stress(:,:,:,:)
DOUBLE PRECISION , ALLOCATABLE :: strainrate(:,:,:,:),wdot(:,:,:)
DOUBLE PRECISION , ALLOCATABLE :: seq(:,:,:), elas_ener(:,:,:)
DOUBLE PRECISION , ALLOCATABLE :: ELAST(:,:,:,:,:)
DOUBLE PRECISION , ALLOCATABLE :: Shear_stress(:,:,:)
DOUBLE PRECISION :: epav(3,3),edotpav(3,3)
DOUBLE PRECISION :: defgradmacroactual(3,3),defgradmacrot(3,3)
DOUBLE PRECISION :: velgradmacro(3,3)
DOUBLE PRECISION :: cg(3,3,3,3),cg66aux(6,6)

! by anand

!integer :: ii,jj,kk,ll,iter,i,j,k,l, K1
INTEGER :: kzz,kyy,kxx,kx,ky,kz,m,n,i,j,k, ig ,k1
DOUBLE PRECISION :: xknorm,det,vm, ZERO, ONE



DOUBLE PRECISION :: aux5(5), press
DOUBLE PRECISION :: auxd(3,3),auxs(3,3)
DOUBLE PRECISION :: aux55(5,5),aux3333(3,3,3,3)
INTEGER :: IC, IS !, JPH
DOUBLE PRECISION :: eel(3,3),cinvg(3,3,3,3)
DOUBLE PRECISION :: intime,fintime,wallclock

! for fftw

    complex ( kind = dp ) , allocatable , dimension ( : , : , : ) :: in_forward !> @var Array input in FFT
    complex ( kind = dp ) , allocatable , dimension ( : , : , : ) :: out_forward !> @var Array output in FFT
    complex ( kind = dp ) , allocatable , dimension ( : , : , : ) :: in_backward !> @var Array input of FFT-1
    complex ( kind = dp ) , allocatable , dimension ( : , : , : ) :: out_backward !> @var Array output for FFT-1


   integer (kind = i4b ) :: plan_forward 
   integer (kind = i4b ) :: plan_backward
   integer :: AllocateStatus

 !  real ,allocatable :: T1(:,:,:,:,:) !,T2(:,:,:,:,:)

 
call cpu_time(intime)

! twinframe

twin_frame(1)=-.3528212
twin_frame(2)= .6111043
twin_frame(3)=-.04548461
twin_frame(4)=0.
twin_frame(5)=0.

! twin frame ends here
print *, plan_forward
print *, 'first'

UR0=0

open(ur0,file='fft.in',status='old')
read(ur0,*)output_prefix
read(ur0,*)npts1,npts2,npts3


!ALLOCATE (data(2*npts1*npts2*npts3))
ALLOCATE (work(3,3,NPTS1,NPTS2,NPTS3),workim(3,3,NPTS1,NPTS2,NPTS3))
ALLOCATE (jgrain(NPTS1,NPTS2,NPTS3),jphase(NPTS1,NPTS2,NPTS3))
ALLOCATE (AG(3,3,NPTS1,NPTS2,NPTS3),CG66(6,6,NPTS1,NPTS2,NPTS3))
ALLOCATE (CRSS(NSYSMX,2,NPTS1,NPTS2,NPTS3),sch(5,NSYSMX,NPTS1,NPTS2,NPTS3))
ALLOCATE (edotp(3,3,NPTS1,NPTS2,NPTS3),ept(3,3,NPTS1,NPTS2,NPTS3))
ALLOCATE (sg(3,3,NPTS1,NPTS2,NPTS3),DEFGRAD(3,3,NPTS1,NPTS2,NPTS3))
ALLOCATE (GAMDOT(NSYSMX,NPTS1,NPTS2,NPTS3))
ALLOCATE (TRIALTAU(NSYSMX,2,NPTS1,NPTS2,NPTS3))
ALLOCATE (VELGRAD(3,3,NPTS1,NPTS2,NPTS3))
ALLOCATE (GACUMGR(NPTS1,NPTS2,NPTS3))
ALLOCATE (RSST(NSYSMX,npts1,npts2,npts3))
!ALLOCATE (STRAINRATE(6,npts1,npts2,npts3)) 
ALLOCATE (STRESS(6,npts1,npts2,npts3)) 
ALLOCATE (SEQ(npts1,npts2,npts3)) 
ALLOCATE (ELAS_ENER(npts1,npts2,npts3)) 
!ALLOCATE (WDOT(npts1,npts2,npts3)) 
ALLOCATE (DISGRAD_CUR(3,3,NPTS1,NPTS2,NPTS3))
ALLOCATE (DISGRAD_PRE(3,3,NPTS1,NPTS2,NPTS3))
!ALLOCATE (DIV_ARRAY(3,NPTS1,NPTS2,NPTS3))
!ALLOCATE (DIV_ARRAY_IM(3,NPTS1,NPTS2,NPTS3))
!ALLOCATE (DIV_REAL(3,NPTS1,NPTS2,NPTS3))
ALLOCATE (ELAST(3,3,NPTS1,NPTS2,NPTS3))
ALLOCATE (SHEAR_STRESS(NPTS1,NPTS2,NPTS3))
!ALLOCATE (T1(3,3,NPTS1+1,NPTS2+1,NPTS3+1))
ALLOCATE (PLAS_SLIP(3,3,NPTS1,NPTS2,NPTS3))


 ALLOCATE (  in_forward(NPTS1, NPTS2, NPTS3),STAT=AllocateStatus)
  IF (AllocateStatus == 0 ) print *, 'Successful allocation of in_forw'
 ALLOCATE ( out_forward(NPTS1, NPTS2, NPTS3),STAT=AllocateStatus)
  IF (AllocateStatus == 0 ) print *, 'Successful allocation of out_forw '
ALLOCATE ( in_backward(NPTS1, NPTS2, NPTS3),STAT=AllocateStatus)
  IF (AllocateStatus == 0 ) print *, 'Successful allocation of in_back'
 ALLOCATE (out_backward(NPTS1, NPTS2, NPTS3),STAT=AllocateStatus)
  IF (AllocateStatus == 0 ) print *, 'Successful allocation of out_back'

 CALL OUTPUT_WRITE(0)

write(21,*) 'IT     ERRE       ERRS       SVM' 

nn(1)=npts1
nn(2)=npts2
nn(3)=npts3
!
nn2(1)=npts1
nn2(2)=npts2
!
prodnn=float(nn(1)*nn(2)*nn(3))
wgt=1./prodnn
!
UR1=1      ! FILECRYS
UR2=2      ! FILETEXT

call vpsc_input 

!if(ILATTICE == 1) call lattice_strains(0)


defgradmacrot(1:3,1:3)=0.
defgradmacro(1:3,1:3)=udot(1:3,1:3)*tdot 


DISGRAD_PRE(1:3,1:3,1:NPTS1,1:NPTS2,1:NPTS3) = 0.d0  !initialzing disgrad_pre to zero in the first step
      EDOTP(1:3,1:3,1:NPTS1,1:NPTS2,1:NPTS3) = 0.d0
        EPT(1:3,1:3,1:NPTS1,1:NPTS2,1:NPTS3) = 0.d0 
	DEFGRAD(1:3,1:3,1:NPTS1,1:NPTS2,1:NPTS3) = 0.d0
  PLAS_SLIP(1:3,1:3,1:NPTS1,1:NPTS2,1:NPTS3) = 0.d0

 !cx
 !cx    ELASTIC INITIAL GUESS FOR THE STRESS
 !cx

do  k=1,npts3
do  j=1,npts2
do  i=1,npts1

   cg66aux(1:6,1:6)=cg66(1:6,1:6,i,j,k)

   call chg_basis(aux6,aux33,cg66aux,cg,3,6)

    do ii=1,3
       do jj=1,3
	       sg(ii,jj,i,j,k)=0.
		   do kk=1,3
		   do ll=1,3
			   sg(ii,jj,i,j,k)=sg(ii,jj,i,j,k)+cg(ii,jj,kk,ll)*defgradmacro(kk,ll)
		   enddo
		   enddo
       enddo
    enddo
     
end do
end do
end do

! stop
!     evm=0.

evm=dvm*tdot

write(55,*)'EVM         EVMP        DVM         DVMP        SVM'  
write(56,'(t2,a,t14,a,t26,a,t38,a,t50,a,t62,a,t74,a,t86,a,t98,a,t110,a,t122,a,t134,a,t146,a,t158,a,t170,a,t182,a, &
     &    t194,a,t206,a,t218,a,t230,a)') 'EVM','E11','E22','E33','EPVM','EP11','EP22','EP33','DVM','D11','D22','D33','DPVM', &
     &    'DP11','DP22','DP33','SVM','S11','S22','S33'

! HERE BEGINS THE START OVER THE LOOPS

DO imicro=1,nsteps

!if (imicro == nsteps .or. mod(imicro,ifreq)==0 .or. imicro .le. 2 )  then
if (imicro .le. 2 .or. mod(imicro,ifreq)==0 .or. imicro .ge. n_inc )  then

    call output_write(1)

end if

    write(*,*) '***************************************************'
    write(*,*) 'STEP = ',imicro


    if(nsteps.ne.1) write(21,*) 'STEP = ',imicro
    write(25,*) 'STEP = ',imicro

     !  if (imicro .ge. 2  ) udot = 0.
	   if (imicro .ge. 2 .and. imicro .le. n_inc) then
                udot = 0.
           elseif (imicro .ge. n_inc)then
             udot = udot_in
           end if 
	
        do k=1,npts3
	do j=1,npts2
	do i=1,npts1
!		defgrad(1:3,1:3,i,j,k)=defgradmacro(1:3,1:3)
        do ii= 1,3   ! this would mean taking last converged defgrad
	 do jj =1,3   ! as intial guess for the current step
    	  defgrad(ii,jj,i,j,k)=defgrad(ii,jj,i,j,k)+udot(ii,jj)*tdot  
	 end do
	 end do
	end do
	end do
	end do


     if(imicro .eq. 1 .or. iupdate .eq. 1) call update_schmid
     if(imicro .eq. 2) call update_schmid
     if(imicro .eq. 2) call update_elasticity

     ddefgradmacroacum(1:3,1:3)=0.
     ddefgradmacro(1:3,1:3)=0. ! AK

      iter=0
      erre=2*error
      errs=2*error 

      div_domain = 5e9

   DO WHILE (iter.lt.itmax.and.(errs.gt.error.or.erre.gt.error))

!  DO WHILE (iter .lt. itmax .and. (Maxval(div_domain) .ge. 200) )


    ! LOOP OVER ITERATIONS IN A GIVEN STEP

	iter=iter+1
    
    write(*,*)'ITER = ',iter
    write(25,*)'ITER = ',iter
    write(*,*) 'DIRECT FFT OF STRESS FIELD'

!a    do ii=1,3
!a    do jj=1,3
	
     
 do i =1,6  
	if (i == 1)then 
	ii=1
	jj=1    
    elseif(i == 2)then 
	ii=1
	jj=2    
	elseif(i == 3)then 
	ii=1
	jj=3    
	elseif(i == 4)then 
	ii=2
	jj=2
	elseif(i == 5)then 
	ii=2
	jj=3    
	elseif(i == 6)then 
	ii=3
	jj=3    
	endif
    
          
		  in_forward(1:NPTS1,1:NPTS2,1:NPTS3) = CMPLX(sg(ii,jj,1:NPTS1,1:NPTS2,1:NPTS3),0.d0)
               
		  
		  if (imicro == 1 .and. iter  == 1 .and. ii == 1 .and. jj == 1 ) then
            call dfftw_plan_dft_3d_(plan_forward,NPTS1,NPTS2,NPTS3,in_forward,out_forward,fftw_forward,fftw_estimate)
            print *, 'successful plan creation for forward'  
          endif    

          call dfftw_execute_(plan_forward)

 		    work(ii,jj,1:NPTS1,1:NPTS2,1:NPTS3) =  DBLE(out_forward(1:NPTS1,1:NPTS2,1:NPTS3))
		  workim(ii,jj,1:NPTS1,1:NPTS2,1:NPTS3) = AIMAG(out_forward(1:NPTS1,1:NPTS2,1:NPTS3))

  end do

! filling up the symmetry parts

        work(2,1,1:NPTS1,1:NPTS2,1:NPTS3)   = work(1,2,1:NPTS1,1:NPTS2,1:NPTS3)
		work(3,1,1:NPTS1,1:NPTS2,1:NPTS3)   = work(1,3,1:NPTS1,1:NPTS2,1:NPTS3)
		work(3,2,1:NPTS1,1:NPTS2,1:NPTS3)   = work(2,3,1:NPTS1,1:NPTS2,1:NPTS3)
        workim(2,1,1:NPTS1,1:NPTS2,1:NPTS3) = workim(1,2,1:NPTS1,1:NPTS2,1:NPTS3)
		workim(3,1,1:NPTS1,1:NPTS2,1:NPTS3) = workim(1,3,1:NPTS1,1:NPTS2,1:NPTS3)
		workim(3,2,1:NPTS1,1:NPTS2,1:NPTS3) = workim(2,3,1:NPTS1,1:NPTS2,1:NPTS3)

! end of filling up the symmetry

           !print *, 'ater outforward'

!
!A      if(npts3.gt.1) then
!A		call fourn(data,nn,3,1)
!A      else
!A		call fourn(data,nn2,2,1)
!A      endif
!
!A      k1=0
!A      do  kzz=1,npts3
!A      do  kyy=1,npts2
!A      do  kxx=1,npts1
!A		  k1=k1+1
!A		  work(ii,jj,kxx,kyy,kzz)=data(k1)
          !     write(876,*)data(k1)
!A		  k1=k1+1
!A		  workim(ii,jj,kxx,kyy,kzz)=data(k1)
           !    write(876,*)data(k1) 
!A      end do
!A      end do
!A      end do
      
!a      end do
!a      end do

 !DIV_ARRAY = 0.
 !DIV_ARRAY_IM = 0.
 !DIV_DOMAIN = 0.

    write(*,*) 'CALCULATING G^pq,ij : SG^ij ...'
!
    do kzz=1,npts3
    do kyy=1,npts2
    do kxx=1,npts1

      if(kxx.le.npts1/2) kx=kxx-1
      if(kxx.gt.npts1/2) kx=kxx-npts1-1

      if(kyy.le.npts2/2) ky=kyy-1
      if(kyy.gt.npts2/2) ky=kyy-npts2-1

      if(kzz.le.npts3/2) kz=kzz-1
      if(kzz.gt.npts3/2) kz=kzz-npts3-1

      xk(1)=kx/(delt(1)*nn(1))
      xk(2)=ky/(delt(2)*nn(2))

      if(npts3.gt.1) then
	      xk(3)=kz/(delt(3)*nn(3))
      else
		  xk(3)=0.
      endif

      xknorm=sqrt(xk(1)**2+xk(2)**2+xk(3)**2)

      if (xknorm.ne.0.) then
		  do i=1,3
		  xk2(i)=xk(i)/(xknorm*xknorm*2.*pi)
      xk(i)=xk(i)/xknorm
		  enddo
      endif

!******* addition of cubic interpolation filter : REf:

!if ( xknorm .ne. 0) then
! sinc = sin(xknorm)/(xknorm)
! filter = sinc**2 - (sin(2*xknorm)/(2*xknorm))
! filter = 3*filter/(xknorm**2)
!   work(1:3,1:3,kxx,kyy,kzz)=work(1:3,1:3,kxx,kyy,kzz)*filter
!   workim(1:3,1:3,kxx,kyy,kzz)=workim(1:3,1:3,kxx,kyy,kzz)*filter
!end if

!******** end of filter addition


! calculating the divergence vector of stress in the 
! fourier space for each frequency

!AKK do ii =1,3
!AKK 	div_array(ii,kxx,kyy,kzz) = (-1.*2.*pi)*workim(ii,1,kxx,kyy,kzz)*xk(1) + &
!AKK	                            (-1.*2.*pi)*workim(ii,2,kxx,kyy,kzz)*xk(2) + &
!AKK								(-1.*2.*pi)*workim(ii,3,kxx,kyy,kzz)*xk(3) 

!AKK	div_array_im(ii,kxx,kyy,kzz) = (2.*pi)*work(ii,1,kxx,kyy,kzz)*xk(1) + &
!AKK	                               (2.*pi)*work(ii,2,kxx,kyy,kzz)*xk(2) + &
!AKK								   (2.*pi)*work(ii,3,kxx,kyy,kzz)*xk(3) 
!AKK end do

! end of divergence vector

!AKK correction to higher frequency terms as propsoed by SUQUET
!AKK   if(kxx.eq.npts1/2+1.or.kyy.eq.npts2/2+1.or.(npts3.gt.1.and.kzz.eq.npts3/2+1)) then
!AKK   g1=-s0
!AKK  else 

!**********************************************************
! SUQUET'S CORRRECTION IS FOR EVEN NUMBER OF GRID POINTS
! SINCE FFTW WE MAKE SURE WE HAVE ODD NUMBER OF GRID POINTS
!**********************************************************


	  do i=1,3
      do k=1,3
      a(i,k)=0.
		  do j=1,3
		  do l=1,3
        a(i,k)=a(i,k)+c0(i,j,k,l)*xk(j)*xk(l)
		  end do
		  end do
	  end do
	  end do

    call minv(a,3,det,minv1,minv2)

      do p=1,3
		  do q=1,3
			  do i=1,3
				  do j=1,3
            g1(p,q,i,j)=-a(p,i)*xk(q)*xk(j)
				  end do
			  end do
		  end do
	  end do
  
!AKK end if  ! end for suquet's correction
	  
      do i=1,3
		  do j=1,3
				  ddefgrad(i,j)=0.
				  ddefgradim(i,j)=0.
				  if(kx.eq.0.and.ky.eq.0.and.kz.eq.0.)  goto 4
				  do k=1,3
					  do l=1,3
						  ddefgrad(i,j)  = ddefgrad(i,j)+g1(i,j,k,l)*work(k,l,kxx,kyy,kzz)
						  ddefgradim(i,j)=ddefgradim(i,j)+g1(i,j,k,l)*workim(k,l,kxx,kyy,kzz)
					  enddo
				  enddo  
				  4 continue !end if 
		  end do
	  end do	

	  do i=1,3
		do j=1,3
			work(i,j,kxx,kyy,kzz)=ddefgrad(i,j)
			workim(i,j,kxx,kyy,kzz)=ddefgradim(i,j)
		enddo
      enddo
	end do ! KXX
	end do ! KYY
	end do ! KZZ

! calculating the norm of the divergence vector by summing individual compoenents of the 
! vector over each voxel in the domain using Plancheral's Theroem

!AKK    do kzz=1,npts3
!AKK    do kyy=1,npts2
!AKK    do kxx=1,npts1
!AKK		div_domain(1:3) = div_domain(1:3) + div_array(1:3,kxx,kyy,kzz)**2+ div_array_im(1:3,kxx,kyy,kzz)**2
!AKK	end do
!AKK	end do
!AKK	end do

!AKK	div_domain(1:3)=sqrt(div_domain(1:3))

!AKK	write(6532,*)div_domain(1:3)

! end of the norm of the divergence


	write(*,*) 'INVERSE FFT TO GET STRAIN FIELD'

    do  m=1,3
    do  n=1,3

        in_backward(1:NPTS1,1:NPTS2,1:NPTS3)=CMPLX(work(m,n,1:NPTS1,1:NPTS2,1:NPTS3),workim(m,n,1:NPTS1,1:NPTS2,1:NPTS3))


	  if (imicro == 1 .and. iter  == 1 .and. m == 1 .and. n == 1 ) then
		  call dfftw_plan_dft_3d_(plan_backward,NPTS1,NPTS2,NPTS3,in_backward,out_backward,fftw_backward,fftw_estimate)
		  print *, 'successful plan creation for backward'  
      endif    

      call dfftw_execute_(plan_backward)

		   do  kzz = 1,NPTS3  
           do  kyy = 1,NPTS2 
           do  kxx = 1,NPTS1 
               defgrad(m,n,kxx,kyy,kzz)=defgrad(m,n,kxx,kyy,kzz)+ddefgradmacro(m,n) & 
			   +  DBLE(out_backward( kxx,kyy,kzz ))/(NPTS1*NPTS2*NPTS3) 
           end do
           end do
           end do

    end do
	end do


!
! SYMMETRIZATION 
!
! EXERCISE CAUTION : UPTIL THIS POINT defgrad was DISPLACMENT GRADIENT 
! HERE WE TAKE THE SYMMETRIC PART OF IT TO CALCULATE TOTAL STRAIN.    
! SO DEFGRAD = TOTAL STRAIN 
	
do kzz=1,npts3
do kyy=1,npts2
do kxx=1,npts1 
	  
      DISGRAD_CUR(1:3,1:3,kxx,kyy,kzz)=defgrad(1:3,1:3,kxx,kyy,kzz)
      
      do ii=1,3
      do jj=1,3
             aux33(ii,jj)=(defgrad(ii,jj,kxx,kyy,kzz)+defgrad(jj,ii,kxx,kyy,kzz))/2.
      enddo
      enddo
!
      do ii=1,3
      do jj=1,3
             defgrad(ii,jj,kxx,kyy,kzz)=aux33(ii,jj)
      enddo
      enddo
!
end do ! kxx
end do ! kyy
end do ! kzz


    write(*,*) 'UPDATE STRESS FIELD'

    call evpal
    call get_smacro

    erre=erre/evm
    errs=errs/svm

    write(*,*) 'STRESS FIELD ERROR =',errs
    write(*,*) 'STRAIN FIELD ERROR =',erre
    write(21,fmt='(i3,4(1x,e10.4),10(1x,F7.4))') iter,erre,errs,svm

 END DO !DO WHILE 



! CALCULATING THE DIVERGENCE VECTOR IN REAL SPACE 
! BY INVERSE FFT OF THE DIV_ARRAY
!akk	write(*,*) 'INVERSE FFT TO GET DIVERGENCE VECTOR'
!akk    do  m=1,3
!akk      k1=0
!akk      do k=1,npts3
!akk      do j=1,npts2
!akk      do i=1,npts1 
!akk		  k1=k1+1
!akk		  data(k1)=div_array(m,i,j,k)
!akk		  k1=k1+1
!akk		  data(k1)=div_array_im(m,i,j,k)
!akk      end do
!akk	  end do
!akk	  end do
!
!akk      if(npts3.gt.1) then
!akk		call fourn(data,nn,3,-1)
!akk      else
!akk		call fourn(data,nn2,2,-1)
!akk      endif
!
!akk      do i=1,2*npts1*npts2*npts3
!akk      data(i)=data(i)/prodnn
!akk      enddo
!
!akk     k1=0
!akk      do k=1,npts3
!akk      do j=1,npts2
!akk      do i=1,npts1	
!akk      k1 = k1+1
!akk    	 div_real(m,i,j,k) = data(k1)
!akk	  k1 = k1+1
!akk      end do
!akk	  end do
!akk	  end do
!akk    end do
! END OF REAL DIVERGENCE VECTOR
! here we calcaulte the velgrad = current disgrad - previous disgrad
! velgrad = velgrad/ tdot
! Velgrad is used for reorientation of the crystals
!akk do kzz = 1,npts3
!akk do kyy = 1,npts2
!akk do kxx = 1,npts1
!akk velgrad(1:3,1:3,kxx,kyy,kzz) = DISGRAD_CUR(1:3,1:3,kxx,kyy,kzz) - DISGRAD_PRE(1:3,1:3,kxx,kyy,kzz)
!akk velgrad(1:3,1:3,kxx,kyy,kzz)=velgrad(1:3,1:3,kxx,kyy,kzz)/tdot
!akk DISGRAD_PRE(1:3,1:3,kxx,kyy,kzz) = DISGRAD_CUR(1:3,1:3,kxx,kyy,kzz)
!akk end do
!akk end do
!akk end do

VELGRAD = DISGRAD_CUR - DISGRAD_PRE
VELGRAD = VELGRAD/TDOT
DISGRAD_PRE = DISGRAD_CUR

do ii=1,3
do jj=1,3
	   defgradmacroactual(ii,jj)=defgradmacro(ii,jj)+ddefgradmacroacum(ii,jj)
	   velgradmacro(ii,jj)=defgradmacroactual(ii,jj)-defgradmacrot(ii,jj)
enddo
enddo

   write(*,*)'defgradmacro(1,1),defgradmacro(2,2),defgradmacro(3,3)'
   write(*,*)defgradmacroactual(1,1),defgradmacroactual(2,2),defgradmacroactual(3,3)
   write(*,*)'defgradmacro(1,1)/defgradmacro(3,3)'
   write(*,*)defgradmacroactual(1,1)/defgradmacroactual(3,3)
   write(*,*)'scauav(1,1),scauav(2,2),scauav(3,3)'
   write(*,*)scauav(1,1),scauav(2,2),scauav(3,3)

   evm=vm(defgradmacroactual)
!      dvm=vm(defgradmacroactual-defgradmacrot)/tdot
   dvm=vm(velgradmacro)/tdot  ! in Ricardo's version it is just vm(velgradmacro)
   defgradmacrot=defgradmacroactual
!
!     INITIAL GUESS OF DEFGRADMACRO AT t+dt ALWAYS ELASTIC     
!
!   do ii=1,3
!   do jj=1,3
!      defgradmacro(ii,jj)=defgradmacro(ii,jj)+udot(ii,jj)*tdot
!   enddo
!   enddo
 
   DEFGRADMACRO = DEFGRADMACRO + UDOT*TDOT
   EPT = EPT + EDOTP*TDOT
   PLAS_SLIP = PLAS_SLIP+edotp*tdot

   epav=0.
   edotpav=0.

   do k=1,npts3
   do j=1,npts2
   do i=1,npts1
     do ii=1,3
     do jj=1,3
        epav(ii,jj)=epav(ii,jj)+ept(ii,jj,i,j,k)*wgt
        edotpav(ii,jj)=edotpav(ii,jj)+edotp(ii,jj,i,j,k)*wgt
     enddo
     enddo
   end do
   end do
   end do 

   evmp=0.
   dvmp=0.

   do ii=1,3
   do jj=1,3
       evmp=evmp+epav(ii,jj)**2
       dvmp=dvmp+edotpav(ii,jj)**2
   enddo
   enddo

   evmp=sqrt((2./3.)*evmp)
   dvmp=sqrt((2./3.)*dvmp)

   write(55,315) evm,evmp,dvm,dvmp,svm 

   if (imicro == 1 ) call twin_reorient 

   IF(IUPDATE.EQ.1) THEN

!c     VELMAX

!c      velmax(1)=dsim(1,1)*delt(1)*(npts1-1)
!c      velmax(2)=dsim(2,2)*delt(2)*(npts2-1)
!c      velmax(3)=dsim(3,3)*delt(3)*(npts3-1)

!ccc     UPDATE ORIENTATIONS

      if(iwtex == 1)call update_orient

!     UPDATE DELT

!c      delt(1)=(delt(1)*(npts1-1)+velmax(1)*tdot)/(npts1-1)
!c      delt(2)=(delt(2)*(npts2-1)+velmax(2)*tdot)/(npts2-1)
!c      if(npts3.gt.1) then
!c      delt(3)=(delt(3)*(npts3-1)+velmax(3)*tdot)/(npts3-1)
!c      endif

      ENDIF                !  IUPDATE ENDIF

!     IF(IUPHARD.EQ.1) CALL HARDEN

   write(56,315) evm,defgradmacroactual(1,1), defgradmacroactual(2,2),defgradmacroactual(3,3), &
     &   evmp,epav(1,1),epav(2,2),epav(3,3),dvm,velgradmacro(1,1),velgradmacro(2,2),velgradmacro(3,3), &
     &   dvmp,edotpav(1,1),edotpav(2,2),edotpav(3,3), svm,scauav(1,1),scauav(2,2),scauav(3,3)

315   format(20(e11.4,1x)) 



!c if (iwtex == 1) then

!     TEX.OUT

!c      open(24,file='tex.out',status='unknown')

!c      zero=0.
!c      one=1.

!c      do i=1,3
!c      do j=1,3
!c      fbar(i,j)=0.
!c      enddo
!c      enddo

!c      do i=1,3
!c      fbar(i,i)=delt(i)
!c      enddo

!c      write(24,*)'TEX.OUT from FFT'
!c      write(24,*)'Formatted to plot with POLE'
!c      write(24,111)((fbar(i,j),j=1,3),i=1,3)
!c      write(24,'(a1,i10)') 'B',npts1*npts2*npts3
!c111   format(9f7.3) 

!c      ig=0


!c   do k=1,npts3
!c   do j=1,npts2
!c   do i=1,npts1

!c      ig=ig+1

 !c     do ii=1,3
 !c     do jj=1,3
 !c     aux33(ii,jj)=ag(jj,ii,i,j,k)
 !c     enddo
 !c     enddo

 !c     call euler(1,ph,th,om,aux33)
!cwc
!cwc     together with the euler angles, 
!cwc     write equivalent strain-rate ...
!cwc
!cw      deqx=0.
!cw      seqx=0.
!cw      do jj=1,5
!cw      deqx=deqx+(dbar(jj)+dtilde(jj,i,j,k))**2
!cw      seqx=seqx+sg(jj,i,j,k)**2
!cw      enddo
!cw      deqx=sqrt(2./3.*deqx)
!cw      seqx=sqrt(3./2.*seqx)
!cwc
!cwc     ... and activities of the 2 first active modes 
!cwc         for each Fourier point
!c
!cww      gmodgr1=gmodgr(1,i,j,k)
!cw      gmodgr2=gmodgr(2,i,j,k)
!c
!cw      gmodgr1=0.
!cw      gmodgr2=0.
!cwc
!cw      write(kfile,172) ph,th,om,one,eeq(i,j,k),
!cw     #  deqx,seqx,zero,ig,jgrain(i,j,k),jphase(i,j,k)
!cw172   format(3f8.2,f6.2,4e11.3,3i8)
!cw
!c      write(24,3305) ph,th,om,one,i,j,k,jgrain(i,j,k),jphase(i,j,k)
!c 3305  format(4(f9.3,1x),5i5)

!c  enddo
!c  enddo
!c  enddo
!c end if ! iwtex == 1 case


	   elas_ener = 0.
	   tot_elas_ener = 0.

 do k=1,npts3
 do j=1,npts2
 do i=1,npts1


	STRESS(1,i,j,k)=Sg(1,1,i,j,k)!auxs(1,1)
	STRESS(2,i,j,k)=Sg(2,2,i,j,k)!auxs(2,2)
	STRESS(3,i,j,k)=Sg(3,3,i,j,k)!auxs(3,3)
	STRESS(4,i,j,k)=Sg(2,3,i,j,k)!auxs(2,3)
	STRESS(5,i,j,k)=Sg(1,3,i,j,k)!auxs(3,1)
	STRESS(6,i,j,k)=Sg(1,2,i,j,k)!auxs(1,2) 

	! calculating elastic strain


      if(igas(jph).eq.0) then
       do ii=1,6
       do jj=1,6
        cg66aux(ii,jj)=cg66(ii,jj,i,j,k)
       enddo
       enddo
      call lu_inverse(cg66aux,6)
      call chg_basis(aux6,aux33,cg66aux,cinvg,3,6)
      do ii=1,3
       do jj=1,3
       eel(ii,jj)=0.
		   do kk=1,3
		   do ll=1,3
		   eel(ii,jj)=eel(ii,jj)+ cinvg(ii,jj,kk,ll)*sg(kk,ll,i,j,k)
		   enddo
		   enddo
       enddo
       enddo
	   elast(1:3,1:3,i,j,k)=eel(1:3,1:3)
       else     !   igas else

       do ii=1,3
       do jj=1,3
	       eel(ii,jj)=-100.
       enddo
       enddo
	   end if

! calculating elastic energy
      elas_ener(i,j,k) = 0. 

	  elas_ener(i,j,k) = Sg(1,1,i,j,k)*elast(1,1,i,j,k)+  &
						 Sg(1,2,i,j,k)*elast(1,2,i,j,k)+  &
	                     Sg(1,3,i,j,k)*elast(1,3,i,j,k)+  &
	                     Sg(2,1,i,j,k)*elast(2,1,i,j,k)+  &
	                     Sg(2,2,i,j,k)*elast(2,2,i,j,k)+  &
	                     Sg(2,3,i,j,k)*elast(2,3,i,j,k)+  &
	                     Sg(3,1,i,j,k)*elast(3,1,i,j,k)+  &
	                     Sg(3,2,i,j,k)*elast(3,2,i,j,k)+  &
						 Sg(3,3,i,j,k)*elast(3,3,i,j,k)
	   elas_ener(i,j,k) = 0.5*elas_ener(i,j,k)

	   tot_elas_ener = tot_elas_ener + elas_ener(i,j,k)

 !	   write(6578,*)elas_ener(i,j,k)


	press = (1./3.)*(Sg(1,1,i,j,k)+Sg(2,2,i,j,k)+Sg(3,3,i,j,k))

    auxs(1:3,1:3) = Sg(1:3,1:3,i,j,k) 
	
	auxs(1,1)  = auxs(1,1) - press
	auxs(2,2)  = auxs(2,2) - press
	auxs(3,3)  = auxs(3,3) - press

 !   auxs(1:3,1:3) = Sg(1:3,1:3,i,j,k) 

	call chg_basis(aux5,auxs,aux55,aux3333,2,5)

      SEQ(i,j,k)= STRESS(1,i,j,k)**2+STRESS(2,i,j,k)**2+ &
     		  	STRESS(3,i,j,k)**2+2*STRESS(4,i,j,k)**2+ &
               2*STRESS(5,i,j,k)**2+2*STRESS(6,i,j,k)**2

	 SEQ(i,j,k)= SQRT(1.5*SEQ(i,j,k))

	DO IS = 1,NSYST(JPH)
		RSST(IS,i,j,k)=0.
		DO IC = 1,5
				RSST(IS,i,j,k)=RSST(IS,i,j,k)+aux5(IC)*sch(IC,IS,i,j,k)
		END DO
	END DO

Shear_stress(i,j,k) = 0.
! shear_stress in the twin frame -- calculation
	DO IC = 1,5
		Shear_stress(i,j,k)=Shear_stress(i,j,k)+aux5(IC)*twin_frame(IC)
	END DO
	
enddo
enddo
enddo

!   tot_elas_ener = sum(elas_ener)
  write(523,*)imicro, tot_elas_ener

! calculating the averages

elas_twin = 0.
elas_parent = 0.
plas_twin = 0.
plas_parent = 0.
stress_twin = 0.
stress_parent = 0.
shear_twin = 0.
shear_parent = 0.
el_ener_twin = 0.
el_ener_parent = 0.

vol_twin =0
vol_parent= 0


do k = 1,npts3
do j = 1,npts2
do i = 1,npts1

if (jgrain (i,j,k) == 1 ) then

vol_parent= vol_parent+1

stress_parent(1) = stress_parent(1)  + sg(1,1,i,j,k)
stress_parent(2) = stress_parent(2)  + sg(2,2,i,j,k)
stress_parent(3) = stress_parent(3)  + sg(3,3,i,j,k)
stress_parent(4) = stress_parent(4)  + sg(3,2,i,j,k)
stress_parent(5) = stress_parent(5)  + sg(1,3,i,j,k)
stress_parent(6) = stress_parent(6)  + sg(1,2,i,j,k)

elas_parent(1) = elas_parent(1)  + elast(1,1,i,j,k)
elas_parent(2) = elas_parent(2)  + elast(2,2,i,j,k)
elas_parent(3) = elas_parent(3)  + elast(3,3,i,j,k)
elas_parent(4) = elas_parent(4)  + elast(3,2,i,j,k)
elas_parent(5) = elas_parent(5)  + elast(1,3,i,j,k)
elas_parent(6) = elas_parent(6)  + elast(1,2,i,j,k)

plas_parent(1) = plas_parent(1)  + ept(1,1,i,j,k)
plas_parent(2) = plas_parent(2)  + ept(2,2,i,j,k)
plas_parent(3) = plas_parent(3)  + ept(3,3,i,j,k)
plas_parent(4) = plas_parent(4)  + ept(3,2,i,j,k)
plas_parent(5) = plas_parent(5)  + ept(1,3,i,j,k)
plas_parent(6) = plas_parent(6)  + ept(1,2,i,j,k)

shear_parent = shear_parent + Shear_stress(i,j,k)
el_ener_parent = el_ener_parent +elas_ener(i,j,k)

elseif (jgrain(i,j,k) ==3) then

vol_twin= vol_twin+1

stress_twin(1) = stress_twin(1)  + sg(1,1,i,j,k)
stress_twin(2) = stress_twin(2)  + sg(2,2,i,j,k)
stress_twin(3) = stress_twin(3)  + sg(3,3,i,j,k)
stress_twin(4) = stress_twin(4)  + sg(3,2,i,j,k)
stress_twin(5) = stress_twin(5)  + sg(1,3,i,j,k)
stress_twin(6) = stress_twin(6)  + sg(1,2,i,j,k)

elas_twin(1) = elas_twin(1)  + elast(1,1,i,j,k)
elas_twin(2) = elas_twin(2)  + elast(2,2,i,j,k)
elas_twin(3) = elas_twin(3)  + elast(3,3,i,j,k)
elas_twin(4) = elas_twin(4)  + elast(3,2,i,j,k)
elas_twin(5) = elas_twin(5)  + elast(1,3,i,j,k)
elas_twin(6) = elas_twin(6)  + elast(1,2,i,j,k)

plas_twin(1) = plas_twin(1)  + ept(1,1,i,j,k)
plas_twin(2) = plas_twin(2)  + ept(2,2,i,j,k)
plas_twin(3) = plas_twin(3)  + ept(3,3,i,j,k)
plas_twin(4) = plas_twin(4)  + ept(3,2,i,j,k)
plas_twin(5) = plas_twin(5)  + ept(1,3,i,j,k)
plas_twin(6) = plas_twin(6)  + ept(1,2,i,j,k)

shear_twin = shear_twin + Shear_stress(i,j,k)
el_ener_twin = el_ener_twin +elas_ener(i,j,k)


end if

end do
end do
end do

elas_twin =       elas_twin/vol_twin
plas_twin =       plas_twin/vol_twin
stress_twin =   stress_twin/vol_twin
shear_twin =     shear_twin/vol_twin
el_ener_twin = el_ener_twin/vol_twin

elas_parent =       elas_parent/vol_parent
plas_parent =       plas_parent/vol_parent
stress_parent =   stress_parent/vol_parent
shear_parent =     shear_parent/vol_parent
el_ener_parent = el_ener_parent/vol_parent

write(510,fmt='(I5,F8.3,1x,F8.6,1x,6(F8.4,1x),6(F9.7,1x),6(F8.5,1x))')imicro,shear_parent,el_ener_parent, & 
                stress_parent,elas_parent(1:6),plas_parent(1:6)
write(511,fmt='(I5,F8.3,1x,F8.6,1x,6(F8.4,1x),6(F9.7,1x),6(F8.5,1x))')imicro,shear_twin,el_ener_twin,stress_twin, &
                elas_twin(1:6),plas_twin(1:6)

!  IF(IWTEX.EQ.1 .and. imicro == nsteps) THEN

!if (imicro == nsteps  .or. mod(IMICRO,IFREQ) == 0 .or. imicro .le. 2) THEN
if (imicro .le. 2 .or. mod(imicro,ifreq)==0 .or. imicro .ge. n_inc )  then


! writing to the paraview file
	
	do i = 1,7 

!C	 stress  - file 92
	if(i == 1) write(92,fmt='(A17)')'SCALARS grain int'
	if(i == 2) write(92,fmt='(A17)')'SCALARS s11 float'
	if(i == 3) write(92,fmt='(A17)')'SCALARS s22 float'
	if(i == 4) write(92,fmt='(A17)')'SCALARS s33 float'
	if(i == 5) write(92,fmt='(A17)')'SCALARS s23 float'
	if(i == 6) write(92,fmt='(A17)')'SCALARS s13 float'
	if(i == 7) write(92,fmt='(A17)')'SCALARS s12 float'

	write(92,fmt='(A20)')'LOOKUP_TABLE default'
	if (i == 1) then 
	do k=1,npts3
	do j=1,npts2
	write(92,fmt='(256(g6.4))')jgrain(1,j,k)!jgrain(1:npts1,j,k)
	end do
	end do
	elseif ( i .gt. 1 .and. i .le.7 ) then
	do k=1,npts3
	do j=1,npts2
	write(92,fmt='(256(F8.2,1X))')STRESS(i-1,1,j,k)!STRESS(i-1,1:npts1,j,k)
	end do
	end do
	end if

!C	 strain  - file 94
	 if(i == 1) write(94,fmt='(A17)')'SCALARS grain int'
	 if(i == 2) write(94,fmt='(A17)')'SCALARS e11 float'
	 if(i == 3) write(94,fmt='(A17)')'SCALARS e22 float'
	 if(i == 4) write(94,fmt='(A17)')'SCALARS e33 float'
	 if(i == 5) write(94,fmt='(A17)')'SCALARS e23 float'
	 if(i == 6) write(94,fmt='(A17)')'SCALARS e13 float'
	 if(i == 7) write(94,fmt='(A17)')'SCALARS e12 float'

	 if(i == 1) write(98,fmt='(A18)')'SCALARS grain int'
	 if(i == 2) write(98,fmt='(A18)')'SCALARS EL11 float'
	 if(i == 3) write(98,fmt='(A18)')'SCALARS EL22 float'
	 if(i == 4) write(98,fmt='(A18)')'SCALARS EL33 float'
	 if(i == 5) write(98,fmt='(A18)')'SCALARS EL23 float'
	 if(i == 6) write(98,fmt='(A18)')'SCALARS EL13 float'
	 if(i == 7) write(98,fmt='(A18)')'SCALARS EL12 float' 

	 if(i == 1) write(99,fmt='(A18)')'SCALARS grain int'
	 if(i == 2) write(99,fmt='(A18)')'SCALARS SL11 float'
	 if(i == 3) write(99,fmt='(A18)')'SCALARS SL22 float'
	 if(i == 4) write(99,fmt='(A18)')'SCALARS SL33 float'
	 if(i == 5) write(99,fmt='(A18)')'SCALARS SL23 float'
	 if(i == 6) write(99,fmt='(A18)')'SCALARS SL13 float'
	 if(i == 7) write(99,fmt='(A18)')'SCALARS SL12 float'
	
	write(94,fmt='(A20)')'LOOKUP_TABLE default'
	write(98,fmt='(A20)')'LOOKUP_TABLE default'
	write(99,fmt='(A20)')'LOOKUP_TABLE default'


	if (i == 1) then 
	do k=1,npts3
	do j=1,npts2
	write(94,fmt='(256(g6.4))')jgrain(1,j,k) !jgrain(1:npts1,j,k)
	write(98,fmt='(256(g6.4))')jgrain(1,j,k) !jgrain(1:npts1,j,k)
	write(99,fmt='(256(g6.4))')jgrain(1,j,k) !jgrain(1:npts1,j,k)
	end do
	end do
	elseif ( i .eq. 2 ) then
	do k=1,npts3
	do j=1,npts2
	write(94,fmt='(256(F8.3,1X))')ept(1,1,1,j,k)!ept(1,1,1:npts1,j,k)
	write(98,fmt='(256(F9.7,1X))')elast(1,1,1,j,k)!elast(1,1,1:npts1,j,k)
	write(99,fmt='(256(F8.3,1X))')plas_slip(1,1,1,j,k)!plas_slip(1,1,1:npts1,j,k)
	end do
	end do
	elseif ( i .eq. 3 ) then
	do k=1,npts3
	do j=1,npts2
	write(94,fmt='(256(F8.3,1X))')ept(2,2,1,j,k) !ept(2,2,1:npts1,j,k)
	write(98,fmt='(256(F9.7,1X))')elast(2,2,1,j,k) !elast(2,2,1:npts1,j,k)
	write(99,fmt='(256(F8.3,1X))')plas_slip(2,2,1,j,k) !plas_slip(2,2,1:npts1,j,k)
	end do
	end do
	elseif ( i .eq. 4 ) then
	do k=1,npts3
	do j=1,npts2
	write(94,fmt='(256(F8.3,1X))')ept(3,3,1,j,k) ! ept(3,3,1:npts1,j,k)
	write(98,fmt='(256(F9.7,1X))')elast(3,3,1,j,k) !elast(3,3,1:npts1,j,k)
	write(99,fmt='(256(F8.3,1X))')plas_slip(3,3,1,j,k) !plas_slip(3,3,1:npts1,j,k)
	end do
	end do
	elseif ( i .eq. 5 ) then
	do k=1,npts3
	do j=1,npts2
	write(94,fmt='(256(F8.3,1X))')ept(2,3,1,j,k) ! ept(2,3,1:npts1,j,k)
	write(98,fmt='(256(F9.7,1X))')elast(2,3,1,j,k) !elast(2,3,1:npts1,j,k)
	write(99,fmt='(256(F8.3,1X))')plas_slip(2,3,1,j,k) !plas_slip(2,3,1:npts1,j,k)
	end do
	end do
	elseif ( i .eq. 6 ) then
	do k=1,npts3
	do j=1,npts2
	write(94,fmt='(256(F8.3,1X))')ept(1,3,1,j,k) ! ept(1,3,1:npts1,j,k)
	write(98,fmt='(256(F9.7,1X))')elast(1,3,1,j,k) !elast(1,3,1:npts1,j,k)
	write(99,fmt='(256(F8.3,1X))')plas_slip(1,3,1,j,k) !plas_slip(1,3,1:npts1,j,k)
	end do
	end do
	elseif ( i .eq. 7 ) then
	do k=1,npts3
	do j=1,npts2
	write(94,fmt='(256(F8.3,1X))')ept(1,2,1,j,k) !ept(1,2,1:npts1,j,k)
	write(98,fmt='(256(F9.7,1X))')elast(1,2,1,j,k) !elast(1,2,1:npts1,j,k)
	write(99,fmt='(256(F8.3,1X))')plas_slip(1,2,1,j,k) !plas_slip(1,2,1:npts1,j,k)
	end do
	end do

	end if
	end do ! end of loop over components


!C for RSS
	write(93,fmt='(A17)')'SCALARS grain int'
	write(93,fmt='(A20)')'LOOKUP_TABLE default'
	do k=1,npts3
	do j=1,npts2
	write(93,fmt='(256(g6.4))')jgrain(1:npts1,j,k)
	end do
	end do

	DO IS = 1,NSYST(JPH)

	if (IS == 1)write(93,fmt='(A20)')'SCALARS RSS1  float'
	if (IS == 2)write(93,fmt='(A20)')'SCALARS RSS2  float'
	if (IS == 3)write(93,fmt='(A20)')'SCALARS RSS3  float'
	if (IS == 4)write(93,fmt='(A20)')'SCALARS RSS4  float'
	if (IS == 5)write(93,fmt='(A20)')'SCALARS RSS5  float'
	if (IS == 6)write(93,fmt='(A20)')'SCALARS RSS6  float'
	if (IS == 7)write(96,fmt='(A20)')'SCALARS RSS7  float'
	if (IS == 8)write(96,fmt='(A20)')'SCALARS RSS8  float'
	if (IS == 9)write(96,fmt='(A20)')'SCALARS RSS9  float'
	if (IS == 10)write(96,fmt='(A20)')'SCALARS RSS10  float'
	if (IS == 11)write(96,fmt='(A20)')'SCALARS RSS11  float'
	if (IS == 12)write(96,fmt='(A20)')'SCALARS RSS12  float'
	if (IS == 13)write(96,fmt='(A20)')'SCALARS RSS13  float'
	if (IS == 14)write(96,fmt='(A20)')'SCALARS RSS14  float'
	if (IS == 15)write(96,fmt='(A20)')'SCALARS RSS15  float'
	if (IS == 16)write(96,fmt='(A20)')'SCALARS RSS16  float'
	if (IS == 17)write(96,fmt='(A20)')'SCALARS RSS17  float'
	if (IS == 18)write(96,fmt='(A20)')'SCALARS RSS18  float'
	if (IS == 19)write(97,fmt='(A20)')'SCALARS RSS19  float'
	if (IS == 20)write(97,fmt='(A20)')'SCALARS RSS20  float'
	if (IS == 21)write(97,fmt='(A20)')'SCALARS RSS21  float'
	if (IS == 22)write(97,fmt='(A20)')'SCALARS RSS22  float'
	if (IS == 23)write(97,fmt='(A20)')'SCALARS RSS23  float'
	if (IS == 24)write(97,fmt='(A20)')'SCALARS RSS24  float'
	if (IS == 25)write(93,fmt='(A20)')'SCALARS RSS25  float'
	if (IS == 26)write(93,fmt='(A20)')'SCALARS RSS26  float'
	if (IS == 27)write(93,fmt='(A20)')'SCALARS RSS27  float'


    if (is .le. 6 ) then 
	write(93,fmt='(A20)')'LOOKUP_TABLE default'
		do k=1,npts3
		do j=1,npts2
		write(93,fmt='(256(F8.2,1X))')abs(RSST(IS,1,j,k)) !abs(RSST(IS,1:npts1,j,k))
		end do
		end do
	elseif (is .ge. 7  .and. is .le. 18) then
	write(96,fmt='(A20)')'LOOKUP_TABLE default'
		do k=1,npts3
		do j=1,npts2
		write(96,fmt='(256(F8.2,1X))')abs(RSST(IS,1,j,k)) !abs(RSST(IS,1:npts1,j,k))
		end do
		end do
	elseif (is .ge. 19  .and. is .le. 24) then
	write(97,fmt='(A20)')'LOOKUP_TABLE default'
		do k=1,npts3
		do j=1,npts2
		write(97,fmt='(256(F8.2,1X))')RSST(IS,1,j,k) !RSST(IS,1:npts1,j,k)
		end do
		end do
    end if

	END DO

!c
!c write the invariants to the file
!c

	write(95,fmt='(A17)')'SCALARS grain int'
	write(95,fmt='(A20)')'LOOKUP_TABLE default'
	do k=1,npts3
	do j=1,npts2
	write(95,fmt='(256(g6.4))')jgrain(1,j,k) !jgrain(1:npts1,j,k)
	end do
	end do

	write(95,fmt='(A19)')'SCALARS EL_EN float'
	write(95,fmt='(A20)')'LOOKUP_TABLE default'
	do k=1,npts3
	do j=1,npts2
	write(95,fmt='(256(F8.6,1X))')ELAS_ENER(1,j,k) !ELAS_ENER(1:npts1,j,k)
	end do
	end do

	write(95,fmt='(A27)')'SCALARS Shear_stress float'
	write(95,fmt='(A20)')'LOOKUP_TABLE default'
	do k=1,npts3
	do j=1,npts2
	write(95,fmt='(256(F8.2,1X))')Shear_stress(1,j,k) !Shear_stress(1:npts1,j,k)
	end do
	end do

ENDIF

END DO ! ON THE STEPS , i.e. IMICRO


 DEALLOCATE (  in_forward, out_forward, in_backward, out_backward )
 DEALLOCATE ( WORK, WORKIM, JGRAIN, JPHASE, AG )
 DEALLOCATE ( CG66, CRSS, SCH, EDOTP, EPT, SG )
 DEALLOCATE ( DEFGRAD, GAMDOT, TRIALTAU, VELGRAD, GACUMGR )
 DEALLOCATE ( RSST, STRESS, SEQ, ELAS_ENER, ELAST )
 DEALLOCATE ( DISGRAD_CUR, DISGRAD_PRE )
 DEALLOCATE ( SHEAR_STRESS, PLAS_SLIP )


!fintime = time()
call cpu_time(fintime)
wallclock = fintime - intime
write(*,*)"wallclock :: ",wallclock

END PROGRAM FFTW_EIGEN_TWIN2a

!
!******************************************************************
!******************************************************************
!

SUBROUTINE fourn(data,nn,ndim,isign)
IMPLICIT NONE

   INTEGER :: isign,ndim,nn(ndim)
      DOUBLE PRECISION :: data(*)
      INTEGER :: i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot
      DOUBLE PRECISION :: tempi,tempr
      DOUBLE PRECISION :: theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do idim=1,ndim
        ntot=ntot*nn(idim)
      end do 
      nprev=1
   do idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do  i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do  i1=i2,i2+ip1-2,2
              do  i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
             end do 
           end do
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
		end do
        ifp1=ip1

2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do i3=1,ifp1,ip1
            do i1=i3,i3+ip1-2,2
              do i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
            end do
          end do
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
        end do
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
 end do
 
return
END
!
!**************************************************************
!
      SUBROUTINE VOIGT(C2,C4,IOPT)

!
! *** TRANSFORMS SECOND ORDER MATRIX C2 INTO FOURTH ORDER TENSOR C4 IF
! *** IOPT=1 AND VICEVERSA IF IOPT=2. IF IOPT=3,TRANSFORMS WITH INV.FACT.
! *** IOPT=4 FOR GO FROM 6x6 TO 3x3x3x3 WITH Aijkl ANTISYMMETRY
!
	IMPLICIT NONE

     double precision :: C2(6,6),C4(3,3,3,3),IJV(6,2),F(6,6)

!      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/

!by anand

INTEGER :: IOPT, I,I1,I2,J,J1,J2

IJV(1:6,1)= (/1,2,3,2,1,1/)
IJV(1:6,2)= (/1,2,3,3,3,2/)

      IF(IOPT.EQ.1) THEN
		  DO I=1,6
			  I1=IJV(I,1)
			  I2=IJV(I,2)
				  DO J=1,6
					  J1=IJV(J,1)
					  J2=IJV(J,2)
					  C4(I1,I2,J1,J2)=C2(I,J)
					  C4(I2,I1,J1,J2)=C2(I,J)
					  C4(I1,I2,J2,J1)=C2(I,J)
					  C4(I2,I1,J2,J1)=C2(I,J)
				  END DO
		  END DO
      ENDIF

      IF(IOPT.EQ.2) THEN
		  DO I=1,6
			  I1=IJV(I,1)
			  I2=IJV(I,2)
			  DO J=1,6
				  J1=IJV(J,1)
				  J2=IJV(J,2)
				  C2(I,J)=C4(I1,I2,J1,J2)
			  END DO
		  END DO
      ENDIF
      
      IF(IOPT.EQ.3) THEN
		  DO I=1,6
			  DO J=1,6
				  F(I,J)=1.
				  IF(I.GT.3) F(I,J)=0.5
				  IF(J.GT.3) F(I,J)=0.5*F(I,J)
			  END DO
		  END DO

		  DO I=1,6
			  I1=IJV(I,1)
			  I2=IJV(I,2)
				  DO J=1,6
					  J1=IJV(J,1)
					  J2=IJV(J,2)
					  C4(I1,I2,J1,J2)=F(I,J)*C2(I,J)
					  C4(I2,I1,J1,J2)=F(I,J)*C2(I,J)
					  C4(I1,I2,J2,J1)=F(I,J)*C2(I,J)
					  C4(I2,I1,J2,J1)=F(I,J)*C2(I,J)
				  END DO
		  END DO
      ENDIF

      IF(IOPT.EQ.4) THEN
		  DO  I=1,6
			  I1=IJV(I,1)
			  I2=IJV(I,2)
			  DO  J=1,6
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
			  END DO
		  END DO
      ENDIF

  RETURN
 END

!
!**************************************************************
!
!c *****************************************************************************
      subroutine euler(iopt,ph,th,tm,a)
!c
!c     CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
!c     MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
!c     A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.
!c     ph,th,om ARE THE EULER ANGLES (in degrees) OF ca REFERRED TO sa.
!c *****************************************************************************
      implicit none
     
	  integer :: iopt
      double precision :: ph,th,tm,pi
	  double precision :: sph,cph,sth,cth,stm,ctm
      double precision :: a(3,3)
      pi=4.*atan(1.d0)


      if(iopt.eq.1) then
        th=acos(a(3,3))
        if(dabs(a(3,3)).ge.0.9999) then
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

!c
!C ********************************************************************
!C     SUBROUTINE VPSC_INPUT      --->      VERSION 31/jan/99
!C
!C     READS CHARACTERISTICS OF THE RUN: # OF PHASES, NAMES OF INPUT FILES,
!C     DEFORMATION TO BE IMPOSED, CONVERGENCE PARAMETERS, ETC.
!C     READS SINGLE CRYSTAL PROPERTIES: DEFORMATION MODES, CRSS, HARDENING
!C     READS CRYSTAL AND MORPHOLOGIC TEXTURES.
!C     INITIALIZES ARRAYS REQUIRED TO RUN VPSC.
!C     OPENS AND CLOSES INPUT FILES.   OPENS OUTPUT FILES.
!C
!C     MODIFIED 21/07/98 by CNT:
!C     INITIALIZATION RELATED TO 'ELEMENTS' IS DONE INSIDE A SINGLE BLOCK.
!C *****************************************************************************
!C

SUBROUTINE VPSC_INPUT

      USE GLOBAL

	  IMPLICIT NONE

      DOUBLE PRECISION :: aux55(5,5),aux66(6,6),aux3333(3,3,3,3)
      DOUBLE PRECISION :: dsimdev(3,3)

      DOUBLE PRECISION :: IJV(6,2)
! by anand

      DOUBLE PRECISION :: DUM1(6),DUM2(3,3),DUM3(6,6),DUM4(3,3,3,3), dummy ! recheck this very carefully
      integer :: iph, i,j,k 


IJV(1:6,1)= (/1,2,3,2,1,1/)
IJV(1:6,2)= (/1,2,3,3,3,2/)


! *********   INITIALIZATION BLOCK   ***************************
!
!      PI=4.*ATAN(1.)
!
!      IREADSET=0      ! used as control in CUBCOMP to open unit

!C     CALCULATES TENSORS OF THE SYMMETRIC BASIS 'B(3,3,6)'
      CALL CHG_BASIS(DUM1,DUM2,DUM3,DUM4,0,6)
!C
!C     SEED FOR RANDOM NUMBER GENERATOR (RAN2) (USED FOR TWINNING AND RX)
      JRAN=-1

      READ(UR0,*) nph

!     THE FOLLOWING REQUIRED FOR SEVERAL ROUTINES WITH 'do iph=iphbot,iphtop'
        iphbot=1
        iphtop=nph
!
      if(nph.gt.nphmx) then
        write(*,'('' number of phases exceeds multiphase dimens !!'')')
        write(*,'('' --> increase parameter NPHMX to'',i5)') nph
        stop
      endif
!c
      READ(UR0,*) ngr
!c
      if(ngr.ne.npts1*npts2*npts3) then
      write(*,*) 'NUMBER OF FPs DECLARED IN FFT.IN = ',ngr
      write(*,*) 'DOES NOT MATCH WITH NPTS1*NPTS2*NPTS3 IN FFT.DIM =', npts1*npts2*npts3
      stop
      endif
!c
!c     RVE DIMENSIONS
!c
      READ(UR0,*) DELT
      DELTVOL3=(DELT(1)*DELT(2)*DELT(3))**(1./3.)
!c      
      READ(UR0,'(a)') prosa
      READ(UR0,'(a)') filetext

!c ***************************************************************************
!c     LOOP OVER PHASES
!c ***************************************************************************

      DO IPH=1,NPH
!c
      READ(UR0,'(a)') prosa
      READ(UR0,*) igas(iph)
      READ(UR0,'(a)') prosa
      READ(UR0,'(a)') filecryspl
!cevp
      READ(UR0,'(a)') filecrysel
!cevp
!c
!c     READS SLIP AND TWINNING MODES FOR THE PHASE
!c
      if(igas(iph).eq.0) then
      OPEN (unit=UR1,file=filecryspl,status='old')
        call data_crystal(iph)
      CLOSE(unit=UR1)
!cevp
      OPEN (unit=UR1,file=filecrysel,status='old')
        call data_crystal_elast(iph)
      CLOSE(unit=UR1)
!cevp
      endif
!c
      ENDDO     ! END OF DATA INPUT LOOP OVER ALL PHASES
!c
!c     READS INITIAL TEXTURE FROM FILETEXT
!cevp  and calculates the local elastic stiffness
!c
      OPEN(unit=UR2,file=filetext,status='old')
        call data_grain
      CLOSE(unit=UR2)
!ccc
!ccc     SEARCH FOR NRSMIN (NEEDED TO GET TAUMAX INSIDE NR SUBROUTINE)
!ccc
!cc      nrsmin=nrs(1,1)
!cc      DO IPH=1,NPH
!cc        do is=1,nsyst(iph)
!cc          if(nrs(is,iph).lt.nrsmin) nrsmin=nrs(is,iph)
!cc        enddo
!cc      ENDDO
!c
!C ****************************************************************************
!C     READ BOUNDARY CONDITIONS ON OVERALL STRESS AND STRAIN-RATE
!C ****************************************************************************
!C
      READ(UR0,'(A)') PROSA
      READ(UR0,'(A)') PROSA
!C
      do i=1,3
        READ(UR0,*) (iudot(i,j),j=1,3)
      enddo
!c
      if(iudot(1,1)+iudot(2,2)+iudot(3,3).eq.2) then
        write(*,*) 'CHECK DIAGONAL BOUNDARY CONDITIONS IUDOT'
        write(*,*) 'CANNOT ENFORCE ONLY TWO DEVIATORIC COMPONENTS'
        stop
      endif
!c
      do i=1,3
      do j=1,3
        if(i.ne.j.and.iudot(i,j)+iudot(j,i).eq.0) then
          write(*,*) 'CHECK OFF-DIAGONAL BOUNDARY CONDITIONS IUDOT'
          stop
        endif
      enddo
      enddo
!c
      READ(UR0,*)
      DO I=1,3
        READ(UR0,*) (UDOT(I,J),J=1,3)
      ENDDO
     
      UDOT_IN = UDOT

!C
!c     SYMMETRIC STRAIN-RATE, ANTISYMMETRIC ROTATION-RATE TENSORS
!c     AND INDICES OF IMPOSED COMPONENTS
!c
      do i=1,3
      do j=1,3
        dsim(i,j)=(udot(i,j)+udot(j,i))/2.
        tomtot(i,j)=(udot(i,j)-udot(j,i))/2.
      enddo
      enddo
!c
      do i=1,3
        idsim(i)=iudot(i,i)
      enddo
!c
      idsim(4)=0
      if(iudot(2,3).eq.1.and.iudot(3,2).eq.1) idsim(4)=1
      idsim(5)=0
      if(iudot(1,3).eq.1.and.iudot(3,1).eq.1) idsim(5)=1
      idsim(6)=0
      if(iudot(1,2).eq.1.and.iudot(2,1).eq.1) idsim(6)=1
!c
!c     WRITES STRAIN RATE DSIM(I,J) IN b-BASIS AS A 5-DIM VECTOR DBAR(K)
!c
      call chg_basis(dbar6,dsim,aux66,aux3333,2,6)

      do i=1,5
      dbar5(i)=dbar6(i)
      enddo

      call chg_basis(dbar5,dsimdev,aux55,aux3333,1,5)

      dvm=0.
      do i=1,3
      do j=1,3
      dvm=dvm+dsimdev(i,j)**2
      enddo
      enddo
      dvm=sqrt(2./3.*dvm)
!cw
!cw      evm=0.
!Cw
      READ(UR0,*)
      READ(UR0,*) iscau(1),iscau(6),iscau(5)
      READ(UR0,*) iscau(2),iscau(4)
      READ(UR0,*) iscau(3)

      do i=1,6
        if(iscau(i)*idsim(i).ne.0.or.iscau(i)+idsim(i).ne.1) then
          WRITE(*,*) ' CHECK BOUNDARY CONDITS ON STRAIN-RATE AND STRESS'
          WRITE(*,'('' IDSIM = '',6I3)') IDSIM
          WRITE(*,'('' ISCAU = '',6I3)') ISCAU
          STOP
        endif
      enddo

      READ(UR0,*)
      READ(UR0,*) scauchy(1,1),scauchy(1,2),scauchy(1,3)
      READ(UR0,*) scauchy(2,2),scauchy(2,3)
      READ(UR0,*) scauchy(3,3)

      scauchy(3,2)=scauchy(2,3)
      scauchy(3,1)=scauchy(1,3)
      scauchy(2,1)=scauchy(1,2)
!c
!c     IF SOME OFF-DIAG COMPS OF SCAUCHY ARE IMPOSED
!c            WRITE THE CORRESP COMP OF SBAR
!c
!cc      if(iscau(4).eq.1) sbar(3)=sqrt(2.d0)*scauchy(2,3)
!cc      if(iscau(5).eq.1) sbar(4)=sqrt(2.d0)*scauchy(1,3)
!cc      if(iscau(6).eq.1) sbar(5)=sqrt(2.d0)*scauchy(1,2)
!c
      READ(UR0,*)
      READ(UR0,*) DUMMY
      READ(UR0,*) ICTRL
!c
      if(ictrl.eq.-1) tdot=dummy
      if(ictrl.eq.0)  tdot=dummy/dvm
      if(ictrl.gt.0)  then
      write(*,*) 'ICTRL>0 NOT IMPLEMENTED'
      stop
      endif
!c 
!c      if(ictrl.gt.0)  eijincr=dummy
!c
!c      if(ictrl.gt.0) then
!c        if(.not.(idsim(ictrl).eq.1.and.
!c     #      dsim(ijv(ictrl,1),ijv(ictrl,2)).ne.0.)) then
!c         write(*,*) 'ICTRL        =',ictrl
!c         write(*,*) 'IDSIM(ICTRL) =',idsim(ictrl)
!c         write(*,*) 'DSIM(ICRTL)  =',dsim(ijv(ictrl,1),ijv(ictrl,2))
!c         write(*,*) 'ICTRL MUST BE -1 TO IMPOSE DSIM*TDOT'
!c         write(*,*) 'OR IT MUST BE 0 TO CONTROL VON MISES INCREMENT'
!c         write(*,*) 'OR IT MUST IDENTIFY A NON-ZERO STRAIN COMPONENT'
!c         stop
!c        endif
!c      endif
!c
      READ(UR0,*)
      READ(UR0,*) NSTEPS
      READ(UR0,*) ERROR
      READ(UR0,*) ITMAX
      READ(UR0,*) IRECOVER

      if(irecover.eq.1) open(50,file='stress.in',status='old',access='sequential',form='unformatted') 
!c
      READ(UR0,*) ISAVE
      READ(UR0,*) IUPDATE
!c
!cw      if(iupdate.eq.1) open(26,file='update.out',status='unknown')
!c
      READ(UR0,*) IUPHARD
      READ(UR0,*) IWTEX
      READ(UR0,*) IWFIELDS
      READ(UR0,*) ILATTICE
      READ(UR0,*) IFREQ
	  READ(UR0,*) N_INC

!cw      READ(UR0,*) FACT2
!cw      FACT2=1.
!cw  
!c
!c     INITIALIZE CRSS AND ACCUM SHEAR FOR GRAINS
!c
!c      do KKK=ngr(iph-1)+1,ngr(iph)
!c
!c        gtotgr(kkk)=0.
!c
!cw      DO IPH=1,NPH
!c
        do kk=1,npts3
        do jj=1,npts2
        do ii=1,npts1

!c
!cevp        gacumgr(ii,jj,kk)=0.
!c
        iph=jphase(ii,jj,kk)
!c
        do i=1,nsyst(iph)
!c
!cw        fact=1.
!cw        if(jphase(ii,jj,kk).eq.2) fact=fact2
!cw
!cw          crss(i,1,ii,jj,kk)=fact*tau(I,1,iph)
!cw          crss(i,2,ii,jj,kk)=fact*tau(I,2,iph)
!cw
!cw          crss(i,1)=fact*tau(I,1,iph)
!cw          crss(i,2)=fact*tau(I,2,iph)

          crss(i,1,ii,jj,kk)=tau(I,1,iph)
          crss(i,2,ii,jj,kk)=tau(I,2,iph)

		TRIALTAU(I,1,II,JJ,KK)=CRSS(I,1,II,JJ,KK)
  		TRIALTAU(I,2,II,JJ,KK)=CRSS(I,2,II,JJ,KK)


!c
        enddo

      enddo
      enddo
      enddo
!c
!cw      ENDDO
!c
      RETURN
      END
!
!************************************************************
!
!C *****************************************************************************
!C     SUBROUTINE DATA_CRYSTAL        --->      VERSION 03/FEB/2000
!C *****************************************************************************

SUBROUTINE DATA_CRYSTAL(IPH)
  USE GLOBAL
  IMPLICIT NONE

DOUBLE PRECISION :: SN(3),SB(3),CDIM(3)
INTEGER :: ISN(12,4),ISB(12,4)
DOUBLE PRECISION :: aux5(5),aux33(3,3),aux55(5,5),aux3333(3,3,3,3)
!REAL :: HSELFX(10),HLATEX(10,10)

! by anand

INTEGER :: iph,kount,nm,jm,iz,js,nsysx,m,im,is,i,j,k

!INTEGER :: modex,nsmx,nrsx,isectwx, MODE(10)
!REAL :: gamd0x,twshx,tau0xf,tau0xb,tau1x,thet0x,thet1x
DOUBLE PRECISION :: SNOR,QNOR,PROD
!REAL :: hselfx(nmodesx),hlatex(nmodesx,nmodesx)
INTEGER :: nrsxaux

!c
      READ(UR1,'(a)') prosa
      READ(UR1,'(a)') icryst(iph)
      READ(UR1,*)     (cdim(i),i=1,3)
      covera=cdim(3)/cdim(1)
      READ(UR1,*)     nmodesx
      READ(UR1,*)     nmodes(iph)
      READ(UR1,*)     (mode(i),i=1,nmodes(iph))
!C
      IF(NMODES(IPH).GT.NMODMX) THEN
        WRITE(*,'('' NMODES IN PHASE'',I3,'' IS'',I3)') IPH,NMODES(IPH)
        WRITE(*,'('' CHANGE PARAMETER NMODMX IN fft.dim'')')
        STOP
      ENDIF
!c
      NTWMOD(iph)=0
      NSYST(iph) =0
      NTWSYS(iph)=0
      KOUNT=1
!c
!c     START READING DEFORMATION MODES FROM FILECRYS
!c
DO nm=1,nmodesx
!c
        READ(UR1,'(a)') prosa
        READ(UR1,*)     modex,nsmx,nrsxaux,gamd0x,twshx,isectwx
        READ(UR1,*)     tau0xf,tau0xb,tau1x,thet0x,thet1x
        READ(UR1,*)     hselfx(nm),(hlatex(nm,jm),jm=1,nmodesx)
!c
!c     SKIPS nsmx LINES IF THE MODE IS NOT IN THE LIST.
        if(modex.ne.mode(kount)) then
          do iz=1,nsmx
            READ(UR1,*)
          enddo
         cycle  !go to 100
        endif
!C
     IF(THET0X.LT.THET1X) THEN
          WRITE(*,'('' INITIAL HARDENING LOWER THAN FINAL HARDENING FOR MODE'',I3,''  IN PHASE'',I3)') KOUNT,IPH
          STOP
      ENDIF
!C
!C     CASE TAU1=0 CORRESPONDS TO LINEAR HARDENING AND IS INDEPENDENT OF TAU0.
!C     AVOID DIVISION BY ZERO
        IF(TAU1X.LE.1.E-6) THEN
          TAU1X=1.E-6
          THET0X=THET1X
        ENDIF
!c
!c     REORDER HARDENING COEFFICIENTS TO ACCOUNT ONLY FOR ACTIVE MODES
        hselfx(kount)=hselfx(nm)
        do i=1,nmodes(iph)
          hlatex(kount,i)=hlatex(nm,mode(i))
        enddo
!c
!c     SYSTEMS GIVEN IN FOUR INDEX NOTATION: HEXAGONALS AND TRIGONALS
!c     SYSTEMS GIVEN IN THREE INDEX NOTATION: CUBIC AND ORTHORHOMBIC
!c
      IF(icryst(iph).EQ.'HEX' .OR. icryst(iph).EQ.'TRI' .OR. icryst(iph).EQ.'hex' .OR. icryst(iph).EQ.'tri') THEN
        DO J=1,NSMX
          READ(UR1,*) (ISN(J,K),K=1,4),(ISB(J,K),K=1,4)
        ENDDO
      ELSE IF(icryst(iph).EQ.'CUB' .OR. icryst(iph).EQ.'ORT'.OR.icryst(iph).EQ.'cub' .OR. icryst(iph).EQ.'ort') THEN
        DO J=1,NSMX
          READ(UR1,*)(ISN(J,K),K=1,3),(ISB(J,K),K=1,3)
        ENDDO
      ELSE
        WRITE(*,'('' CANNOT IDENTIFY THE CRYSTAL SYMMETRY OF PHASE '',I3)') IPH
        STOP
      ENDIF
!C
      NSM(kount,iph)=nsmx
      IF(TWSHX.NE.0) NTWMOD(iph)=NTWMOD(iph)+1
!C
      IF(NTWMOD(IPH).GT.NTWMMX) THEN
        WRITE(*,'('' NTWMOD IN PHASE'',I3,'' IS'',I3)') IPH,NTWMOD(IPH)
        WRITE(*,'('' CHANGE PARAMETER NTWMMX IN VPSC.DIM'')')
        STOP
      ENDIF
!C
  DO  JS=1,NSM(kount,iph)
      NSYST(iph)=NSYST(iph)+1
      NSYSX=NSYST(iph)
      IF(TWSHX.NE.0) NTWSYS(iph)=NTWSYS(iph)+1
      IF(NSYST(IPH).GT.NSYSMX) THEN
        WRITE(*,'('' NSYST IN PHASE'',I3,'' IS'',I3)') IPH,NSYST(IPH)
        WRITE(*,'('' CHANGE PARAMETER NSYSMX IN VPSC.DIM'')')
        STOP
      ENDIF
!C
!C   DEFINES RATE SENSITIVITY AND CRSS FOR EACH SYSTEM IN THE MODE
!C
      GAMD0(NSYSX,iph) =GAMD0X
      NRS(NSYSX,iph)   =NRSXaux
      TWSH(NSYSX,iph)  =TWSHX
      TAU(NSYSX,1,iph) =TAU0XF
      TAU(NSYSX,2,iph) =TAU0XB
      TAU(NSYSX,3,iph) =TAU1X
      THET(NSYSX,1,iph)=THET0X
      THET(NSYSX,2,iph)=THET1X
!c
      isectw(NSYSX,iph)=isectwx
!C
      IF(icryst(iph).EQ.'HEX' .OR. icryst(iph).EQ.'TRI' .OR. icryst(iph).EQ.'hex' .OR. icryst(iph).EQ.'tri') THEN
        SN(1)= ISN(JS,1)
        SN(2)=(ISN(JS,1)+2.*ISN(JS,2))/SQRT(3.)
        SN(3)= ISN(JS,4)/COVERA
        SB(1)= 3./2.*ISB(JS,1)
        SB(2)=(ISB(JS,1)/2.+ISB(JS,2))*SQRT(3.)
        SB(3)= ISB(JS,4)*COVERA
      ELSE IF(icryst(iph).EQ.'CUB' .OR. icryst(iph).EQ.'ORT'.OR.icryst(iph).EQ.'cub' .OR. icryst(iph).EQ.'ort') THEN
        DO M=1,3
          SN(M)=ISN(JS,M)/CDIM(M)
          SB(M)=ISB(JS,M)*CDIM(M)
        ENDDO
      ENDIF
!C
!C *** NORMALIZES SYSTEM VECTORS AND CHECKS NORMALITY
!C
      SNOR=SQRT(SN(1)*SN(1)+SN(2)*SN(2)+SN(3)*SN(3))
      QNOR=SQRT(SB(1)*SB(1)+SB(2)*SB(2)+SB(3)*SB(3))
      PROD=0.
      DO J=1,3
        DNCA(J,NSYSX,iph)=SN(J)/SNOR
        DBCA(J,NSYSX,iph)=SB(J)/QNOR
        IF(ABS(DNCA(J,NSYSX,iph)).LT.1.E-03) DNCA(J,NSYSX,iph)=0.
        IF(ABS(DBCA(J,NSYSX,iph)).LT.1.E-03) DBCA(J,NSYSX,iph)=0.
        PROD=PROD+DNCA(J,NSYSX,IPH)*DBCA(J,NSYSX,IPH)
      ENDDO
      IF(PROD.GE.1.E-3) THEN
        WRITE(*,'(''SYSTEM'',I4,''  IN MODE'',I4,'' IN PHASE'',I4,''  IS NOT ORTHOGONAL !!'')') JS,NM,IPH
        STOP
      ENDIF
!C
!C   DEFINE SCHMID VECTOR IN CRYSTAL AXES FOR EACH SYSTEM
!C
      DO I=1,3
      DO J=1,3
        AUX33(i,j)=(DNCA(i,NSYSX,IPH)*DBCA(j,NSYSX,IPH)+DNCA(j,NSYSX,IPH)*DBCA(i,NSYSX,IPH))/2.
      ENDDO
      ENDDO
!c
      call chg_basis(aux5,aux33,aux55,aux3333,2,5)
!c
      DO I=1,5
        SCHCA(I,NSYSX,IPH)=AUX5(I)
      ENDDO
   END DO   ! END OF LOOP OVER A GIVEN DEFORMATION MODE

      kount=kount+1

END DO    ! END OF LOOP OVER ALL MODES IN A GIVEN PHASE

!C     INITIALIZE SELF & LATENT HARDENING COEFS FOR EACH SYSTEM OF THE PHASE.
!C     ABSOLUTE UNITS ARE ACCOUNTED FOR BY MODULATING FACTOR IN HARDENING LAW.

      I=0
      DO IM=1,NMODES(iph)
      DO IS=1,NSM(IM,iph)
        I=I+1
        J=0
        DO JM=1,NMODES(iph)
        DO JS=1,NSM(JM,iph)
          J=J+1
          HARD(I,J,IPH)=HLATEX(IM,JM)
        ENDDO
        ENDDO
        HARD(I,I,IPH)=HSELFX(IM)
      ENDDO
      ENDDO

!cc      do i=1,2
!cc      write(*,*) (hard(i,j,iph),j=1,2)
!cc      enddo
!cc      pause

!C     LATENT HARDENING OF SLIP AND TWINNING BY TWINNING IS BASED ON THE
!C     RELATIVE DIRECTIONS OF THE SHEAR DIRECTION AND THE TWIN PLANE
!C     THIS APPROACH IS STILL BEING TESTED (30/4/99)

!C     NSLSYS=NSYST(IPH)-NTWSYS(IPH)
!C     DO IS=1,NSYST(IPH)
!C       DO JT=NSLSYS+1,NSYST(IPH)
!C         IF(IS.NE.JT) THEN
!C           COSA=DBCA(1,IS,IPH)*DNCA(1,JT,IPH)+
!C    #           DBCA(2,IS,IPH)*DNCA(2,JT,IPH)+
!C    #           DBCA(3,IS,IPH)*DNCA(3,JT,IPH)
!C           COSA=ABS(COSA)
!C           HARD(IS,JT,IPH)=HARD(IS,JT,IPH)*(0.5+1.0*COSA)
!C         ENDIF
!C       ENDDO
!C     ENDDO
!C
!C     WRITE(10,'(''  HARDENING MATRIX FOR PHASE'',I3)') IPH
!C     DO I=1,NSYST(IPH)
!C       WRITE(10,'(24F5.1)') (HARD(I,J,IPH),J=1,NSYST(IPH))
!C     ENDDO

!C     VERIFICATION OF TWINNING DATA TO BE SURE PROGRAM WILL RUN PROPERLY

      IF (NMODES(IPH) .GT. 1) THEN
        DO I=2,NSYST(IPH)
          IF(TWSH(I,IPH).EQ.0. .AND. TWSH(I-1,IPH).NE.0.) THEN
            WRITE(*,*) ' WARNING! THE TWINNING MODES MUST FOLLOW THE'
            WRITE(*,*) ' SLIP MODES   -->   REORDER CRYSTAL FILE'
            STOP
          ENDIF
        ENDDO
      ENDIF
!C
      RETURN
      END
!C
!C
!C
!C
!C *****************************************************************************
!C     SUBROUTINE DATA_GRAIN        --->      VERSION 31/mar/99
!C *****************************************************************************

SUBROUTINE DATA_GRAIN

   USE GLOBAL
	IMPLICIT NONE

     double precision :: AA(3,3)

!cx      dimension caux3333(3,3,3,3),caux66(6,6),c066(6,6),s066(6,6)
     double precision :: caux3333(3,3,3,3),caux66(6,6),s066(6,6)
!cx
     double precision :: aux6(6),aux33(3,3),dum
!by anand

integer :: kkk,jgr,i1,j1,k1,l1,i2,j2,k2,l2,i,j,k
!cx
      do i=1,6
      do j=1,6
       c066(i,j)=0.
      enddo
      enddo

 DO kkk=1,ngr

      READ(UR2,*) ph,th,om,ii,jj,kk,jgr,jph

      jgrain(ii,jj,kk)=jgr
      jphase(ii,jj,kk)=jph

!C     CALCULATES THE TRANSFORMATION MATRIX AA WHICH TRANSFORMS FROM
!C     SAMPLE TO CRYSTAL. STORES AG, WHICH TRANSFORMS FROM CRYSTAL TO SAMPLE.

        CALL EULER(2,ph*pi/180,th*pi/180,om*pi/180,aa)

        DO J=1,3
        DO K=1,3
          AG(J,K,ii,jj,kk)=AA(K,J)
        ENDDO
        ENDDO

		!write(65,fmt='(9(f8.4,1x))')aa(1,1:3),aa(2,1:3),aa(3,1:3)
		!write(65,*)aa






!cw
!cw        DO J=1,3
!cw        DO K=1,3
!cw          AG(J,K,jgrain(ii,jj,kk))=AA(K,J)
!cw        ENDDO
!cw        ENDDO
!cw

      do i1=1,3
	  do j1=1,3
	  do k1=1,3
	  do l1=1,3
		  dum = 0. 
			 do i2=1,3
			 do j2=1,3
			 do k2=1,3
			 do l2=1,3
				  dum=dum+aa(i2,i1)*aa(j2,j1)*aa(k2,k1) &
							  *aa(l2,l1)*cc(i2,j2,k2,l2)
			! if (i1 == 1 .and. j1==1 .and. k1 == 2 .and. l1 ==3)write(567,*)dum
			 end do	
			 end do	
			 end do	
			 end do	
			 caux3333(i1,j1,k1,l1)=dum
	  end do
	  end do	
	  end do	
	  end do	
!
     !do i1=1,3;do j1=1,3; do k1=1,3;do l1=1,3
     ! write(321,fmt='(4I3,3g12.5)')i1,j1,k1,l1,caux3333(i1,j1,k1,l1),cc(i1,j1,k1,l1)
     !end do; end do;end do;end do
     
      call chg_basis(aux6,aux33,caux66,caux3333,4,6)
!c

        !write(66,fmt='(6g14.5)')caux66(1,1:6)
        !write(66,fmt='(6g14.5)')caux66(2,1:6)
        !write(66,fmt='(6g14.5)')caux66(3,1:6)
        !write(66,fmt='(6g14.5)')caux66(4,1:6)
        !write(66,fmt='(6g14.5)')caux66(5,1:6)
        !write(66,fmt='(6g14.5)')caux66(6,1:6)

     !stop
     do i=1,6 
     do j=1,6
      cg66(i,j,ii,jj,kk)=caux66(i,j)
      c066(i,j)=c066(i,j)+caux66(i,j)*wgt
      enddo
      enddo

 ENDDO

       print"('initialize_etc(): Final c066 : ',6(6(1x,g12.4)/))"  &
            , ((c066(i,j),i=1,6),j=1,6)

     !stop
 
      s066=c066
      call lu_inverse(s066,6)      

      call chg_basis(aux6,aux33,c066,c0,3,6)
      call chg_basis(aux6,aux33,s066,s0,3,6)

    RETURN
END

!C ************************************************************************

!C ************************************************************************
!C     SUBROUTINE CHG_BASIS    --->   VERSION 19/JUL/01
!C
!C     (modif. 06/FEB/98 - same convention as SELFPOLY - C.N.T.)
!C     (modif. 16/JUN/99 - same convention as Maudlin  - C.N.T.)
!C     (modif. 10/MAY/01 - KDIM version - R.L.)
!C
!C     KDIM=5 or 6, FOR DEVIATORIC or DEV+HYDROST TENSORS, RESPECTIVELY.
!C     IOPT=0: DEFINES A BASIS OF 6 SECOND ORDER TENSORS B(N).
!C     IOPT=1: CALCULATES SECOND ORDER TENSOR 'C2' AS AN EXPANSION IN TERMS
!C             OF VECTOR COMPONENTS CE2(KDIM) AND THE BASIS TENSORS B(KDIM).
!C     IOPT=2: CALCULATES COMPONENTS OF C2 AS A VECTOR CE2(KDIM).
!C     IOPT=3: CALCULATES FOURTH ORDER TENSOR 'C4' AS AN EXPANSION IN TERMS
!C             OF MATRIX COMPONENTS CE4(K,K) AND THE BASIS TENSORS B(KDIM).
!C     IOPT=4: CALCULATES MATRIX COMPONENTS CE4(K,K) OF TENSOR 'C4'.
!C **************************************************************************

SUBROUTINE CHG_BASIS(CE2,C2,CE4,C4,IOPT,KDIM)

IMPLICIT NONE

DOUBLE PRECISION :: RSQ2, RSQ3, RSQ6
!c      PARAMETER (SQR2=1.41421356237309   )
!      PARAMETER (RSQ2=0.70710678118654744)
!      PARAMETER (RSQ3=0.57735026918962584)
!      PARAMETER (RSQ6=0.40824829046386304)

INTEGER :: KDIM ! added by anand
DOUBLE PRECISION :: CE2(KDIM),C2(3,3),CE4(KDIM,KDIM),C4(3,3,3,3)
DOUBLE PRECISION :: B(3,3,6)

! by anand

integer :: IOPT,I,J,K,N,L,M,IS
!REAL ::
!C     DIMENSION B(3,3,6)
!C     DATA B /RSQ6,0,   0,   0,   RSQ6,0,   0,   0,  -2*RSQ6,
!C    #        RSQ2,0,   0,   0,  -RSQ2,0,   0,   0,   0,
!C    #        0,   0,   0,   0,   0,   RSQ2,0,   RSQ2,0,
!C    #        0,   0,   RSQ2,0,   0,   0,   RSQ2,0,   0,
!C    #        0,   RSQ2,0,   RSQ2,0,   0,   0,   0,   0,
!C    #        RSQ3,0,   0,   0,   RSQ3,0,   0,   0,   RSQ3/

!      COMMON/BASIS/ B(3,3,6)

RSQ2=0.70710678118654744
RSQ3=0.57735026918962584
RSQ6=0.40824829046386304

!AKK IF(IOPT.EQ.0) THEN
! anand changes this so that when compilers like gfortran do not 
! necessary remember from previous step
! B-basis is calculated again.

IF(IOPT.GE.0) THEN
!C *** CALCULATES BASIS TENSORS B(N)

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

!C *** CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
      IF(IOPT.EQ.1) THEN
        DO  I=1,3
	        DO J=1,3
				C2(I,J)=0.0
				DO N=1,KDIM
				   C2(I,J)=C2(I,J)+CE2(N)*B(I,J,N)
				END DO
			END DO
		END DO
      ENDIF

!C *** CALCULATES KDIMx1 b-COMPONENTS VECTOR FROM SECOND ORDER TENSOR.
      IF(IOPT.EQ.2) THEN
        DO N=1,KDIM
        CE2(N)=0.0
			DO I=1,3
				DO J=1,3
				  CE2(N)=CE2(N)+C2(I,J)*B(I,J,N)
				END DO
			END DO
	    END DO
      ENDIF

!C *** CALCULATES FOURTH ORDER TENSOR FROM b-COMPONENTS MATRIX.
      IF(IOPT.EQ.3) THEN
        DO I=1,3
			DO J=1,3
				DO K=1,3
					DO L=1,3
					C4(I,J,K,L)=0.0
						DO  N=1,KDIM
							DO  M=1,KDIM
							  C4(I,J,K,L)=C4(I,J,K,L)+CE4(N,M)*B(I,J,N)*B(K,L,M)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
      ENDIF

!C *** CALCULATES KDIMxKDIM b-COMPONENTS MATRIX FROM FOURTH ORDER TENSOR.
      IF(IOPT.EQ.4) THEN
        DO N=1,KDIM
			DO M=1,KDIM
			CE4(N,M)=0.0
				DO I=1,3
					DO J=1,3
						DO K=1,3
							DO L=1,3
							  CE4(N,M)=CE4(N,M)+C4(I,J,K,L)*B(I,J,N)*B(K,L,M)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
      ENDIF

      RETURN
  END
!
!******************************************************************************
!
!C **************************************************************************
!C     SUBROUTINE UPDATE_SCHMID
!C
!C     ROTATES SCHMID TENSORS OF EACH GRAIN FROM CRYSTAL TO SAMPLE AXES
!C **************************************************************************
 SUBROUTINE UPDATE_SCHMID

   USE GLOBAL
   IMPLICIT NONE
 
  DOUBLE PRECISION :: aux5(5),aux33(3,3),aux55(5,5),aux3333(3,3,3,3)
  DOUBLE PRECISION :: aux33r(3,3)
! by anand

integer :: is,i1,j1,i,j,k ! jph,


!write(557,*)'23'
DO kk=1,npts3
DO jj=1,npts2
DO ii=1,npts1

!write(345,fmt='(9(f9.5,1x))')ag(1,1:3,ii,jj,kk),ag(2,1:3,ii,jj,kk),ag(3,1:3,ii,jj,kk)

   jph=jphase(ii,jj,kk)

   DO is=1,nsyst(jph)

      do j=1,5
        aux5(j)=schca(j,is,jph)
      enddo

      call chg_basis(aux5,aux33,aux55,aux3333,1,5)

      do i=1,3
		  do j=1,3
			aux33r(i,j)=0.
			  do i1=1,3
				  do j1=1,3
					aux33r(i,j)=aux33r(i,j)+ag(i,i1,ii,jj,kk)*ag(j,j1,ii,jj,kk)*aux33(i1,j1)
				  end do	
			  end do	
		  end do	
	  end do	

      call chg_basis(aux5,aux33r,aux55,aux3333,2,5)

      do j=1,5
        sch(j,is,ii,jj,kk)=aux5(j)
      enddo
!write(345,fmt='(9(f9.5,1x))')sch(1:5,is,ii,jj,kk)

   END DO
ENDDO
ENDDO
ENDDO

return
end
!
!************************************************************************************
!C **************************************************************************
!C     SUBROUTINE UPDATE_ELASTICITY
!C
!C     UPDATE ELASTICITY TENSORS BASED ON NEW CRYSTAL AXES
!C **************************************************************************
 SUBROUTINE UPDATE_ELASTICITY

   USE GLOBAL
   IMPLICIT NONE
 
  DOUBLE PRECISION :: aux5(5),aux33(3,3),aux55(5,5),caux3333(3,3,3,3)
  DOUBLE PRECISION :: aux33r(3,3),cg66aux(6,6),caux66(6,6), cg(3,3,3,3)
  DOUBLE PRECISION :: dum, aa(3,3), aux6(3,3)
! by anand

integer :: is,i1,j1,k1,l1,i2,j2,k2,l2,i,j,k ! jph,k1,


!write(557,*)'23'
DO k=1,npts3
DO j=1,npts2
DO i=1,npts1

!write(345,fmt='(9(f9.5,1x))')ag(1,1:3,ii,jj,kk),ag(2,1:3,ii,jj,kk),ag(3,1:3,ii,jj,kk)


aa(1:3,1:3) = ag(1:3,1:3,i,j,k)

   jph=jphase(i,j,k)

      if(igas(jph).eq.0) then

! change the cg from 66 to 3333

!       do ii=1,6
!       do jj=1,6
!        cg66aux(ii,jj)=cg66(ii,jj,i,j,k)
!       enddo
!       enddo

!      call chg_basis(aux6,aux33,cg66aux,cg,3,6)

! rotate the cc to the new crystal frame

      do i1=1,3
		  do j1=1,3
			  do k1=1,3
				  do l1=1,3
				  DUM =0. 
					  do i2=1,3
						  do j2=1,3
							  do k2=1,3
								  do l2=1,3
									  dum=dum+aa(i2,i1)*aa(j2,j1)*aa(k2,k1)*aa(l2,l1)*cc(i2,j2,k2,l2)
								  end do	
							  end do	
						  end do	
					  end do	
					 caux3333(i1,j1,k1,l1)=dum
				  end do
			  end do	
		  end do	
	  end do	
!
      call chg_basis(aux6,aux33,caux66,caux3333,4,6)
!c
      do ii=1,6
		  do jj=1,6
			  cg66(ii,jj,i,j,k)=caux66(ii,jj)
		!	  c066(ii,jj)=c066(ii,jj)+caux66(ii,jj)*wgt
		  enddo
      enddo
   
   end if

ENDDO
ENDDO
ENDDO

return
end
!
!************************************************************************************
!
! *************************************************************************
      SUBROUTINE LU_INVERSE(A,N)
!
!   INVERTS A MATRIX USING LU DECOMPOSITION
!     
      IMPLICIT NONE
      INTEGER :: I,J,N, ISINGULAR
      DOUBLE PRECISION :: A(N,N),Y(N,N),D
      INTEGER :: INDX(N)
!
      DO I=1,N
        DO J=1,N
          Y(I,J)=0.
        ENDDO
        Y(I,I)=1.
      ENDDO
!
      CALL LUDCMP(A,N,N,INDX,D,ISINGULAR)
!
      DO J=1,N
        CALL LUBKSB(A,N,N,INDX,Y(1,J))
      ENDDO
!
      DO I=1,N
      DO J=1,N
       A(I,J)=Y(I,J)
      ENDDO
      ENDDO
!
      RETURN
      END
!
! *************************************************************************
!
    SUBROUTINE LU_EQSYSTEM(A,B,N,ISINGULAR)
!
!     SOLVES A*X=B USING LU DECOMPOSITION 

      IMPLICIT NONE

      INTEGER :: I,J, ISINGULAR, N , INDX(N)
      DOUBLE PRECISION :: A(N,N),B(N),D

      CALL LUDCMP(A,N,N,INDX,D,ISINGULAR)

      IF(ISINGULAR.EQ.1) RETURN

      CALL LUBKSB(A,N,N,INDX,B)

      RETURN
      END

!C *****************************************************************************
!c
      SUBROUTINE ludcmp(a,n,np,indx,d,isingular)
     
      IMPLICIT NONE

      INTEGER n,np,indx(n),NMAX
!c      REAL d,a(np,np),TINY
!c      PARAMETER (NMAX=500,TINY=1.0e-20)
      DOUBLE PRECISION :: d,a(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,imax,j,k,isingular
      DOUBLE PRECISION :: aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
!c
!c        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
!c
        if(aamax.eq.0.) then
        isingular=1
        return
        endif
!c
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.

        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
!c
!c        if(a(j,j).eq.0.) a(j,j)=TINY
!c
        if(a(j,j).eq.0.) then
        isingular=1
        return
        endif
!c
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
!c
      isingular=0
!c
      return
      END
!c
!c *****************************************************************************
!c
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      DOUBLE PRECISION :: a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION :: sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
!c
!
!***********************************************************************
!
!C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C
!C     FUNCTION TMISMATCH   ---->   VERSION OF 27/DEC/98
!C
!C     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO NROWSxNCOLS MATRICES
!C     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE
!C     OF BOTH DATA.
!C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C
      FUNCTION tmismatch (v1,v2,nrows,ncols)
!C
!COLD      DIMENSION v1(nrows,ncols),v2(nrows,ncols)
!COLD      DIMENSION v_dif(6,6),v_ave(6,6)
      DOUBLE PRECISION :: v1(36),v2(36),tnorm
      DOUBLE PRECISION :: v_dif(36),v_ave(36)
      INTEGER :: I, nrows, ncols

!C    
      do i=1,nrows*ncols
        v_dif(i)=v1(i)-v2(i)
        v_ave(i)=0.5d0*(v1(i)+v2(i))
      enddo
      tmismatch=tnorm(v_dif,nrows,ncols)/tnorm(v_ave,nrows,ncols)
!
      RETURN
      END

!
!C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C
!C     FUNCTION TNORM   ---->   VERSION OF 27/DEC/98
!C
!C     CALCULATES THE NORM OF A NROWSxNCOLS-MATRIX (NROWS,NCOLS =< 6)
!C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!C
      FUNCTION tnorm(v,nrows,ncols)
!C
!COLD      DIMENSION v(nrows,ncols)
      DOUBLE PRECISION :: v(36), tnorm
	  INTEGER :: I, nrows, ncols
!C
      tnorm=0.d0
      do i=1,nrows*ncols
        tnorm=tnorm+v(i)*v(i)
      enddo
      tnorm=sqrt(tnorm)
!C
      RETURN
      END
!
!******************************************************************************
!
!      SUBROUTINE MINV (A,N,D,L,M)
SUBROUTINE MINV (A,N,D,L,M)

  implicit none
  DOUBLE PRECISION :: d , biga , hold 
  integer :: n , i , ij , ik , iz , j , ji , jk , jp , jq
  integer :: jr , k , ki , kj , kk , nk
! logical :: mprint

  !     DIMENSION A(*),L(*),M(*)
  DOUBLE PRECISION :: A(*),L(*),M(*)

  !     SEARCH FOR LARGEST ELEMENT

!  if(mprint) print*,'MINV: A=',(a(i),i=1,n*n)
!  if(mprint) print*,'MINV: N= ',n
  !      D=1.d0
  D=1.0
  NK=-N
  DO 180 K=1,N
     NK=NK+N
     L(K)=K
     M(K)=K
     KK=NK+K
     BIGA=A(KK)
     DO J=K,N
        IZ=N*(J-1)
        DO 20 I=K,N
           IJ=IZ+I
           !               IF (ABS(BIGA)-ABS(A(IJ))) 10,20,20
           ! 10            BIGA=A(IJ)
           !               L(K)=I
           !               M(K)=J
           IF ( ABS(BIGA)-ABS(A(IJ)) .lt. 0. ) then
10            BIGA=A(IJ)
              L(K)=I
              M(K)=J
           end if
20      end do
     end do

     !     INTERCHANGE ROWS

     J=L(K)
     !         IF (J-K) 50,50,30
     IF ( (J-K) .gt. 0 ) then
30      KI=K-N
        DO 40 I=1,N
           KI=KI+N
           HOLD=-A(KI)
           JI=KI-K+J
           A(KI)=A(JI)
           A(JI)=HOLD
40      end do
     end if
     !     INTERCHANGE COLUMNS

50   I=M(K)
     !         IF (I-K) 80,80,60
     IF ( (I-K) .gt. 0 ) then
60      JP=N*(I-1)
        DO 70 J=1,N
           JK=NK+J
           JI=JP+J
           HOLD=-A(JK)
           A(JK)=A(JI)
           A(JI)=HOLD
70      end do
     end if

     !     DIVIDE COLUMN BY MINUS PIVOT (BIGA)

80   IF (ABS(BIGA).LT.1.e-10) THEN
90      D = 0.0
        RETURN
     ENDIF
100  DO 120 I=1,N
        !            print*,'MINV:i,k,nk,ik ',i,k,nk,ik
        !            IF (I-K) 110,120,110
        IF ( (I-K) .ne. 0 ) then
110        IK = NK + I
           !     print*,'MINV-after110: i,k,nk,ik ',i,k,nk,ik
           A(IK) = A(IK) / (-1.*BIGA)
        end if
120  end do

     !     REDUCE MATRIX

     DO I=1,N
        IK=NK+I
        HOLD=A(IK)
        IJ=I-N
        DO J=1,N
           IJ=IJ+N
           !               IF (I-K) 130,150,130
           IF ( (I-K) .ne. 0 ) then
              ! 130              IF (J-K) 140,150,140
              IF ( (J-K) .ne. 0 ) then
140              KJ = IJ-I+K
                 A(IJ) = HOLD*A(KJ)+A(IJ)
              end if
           end if
150     end do
     end do

     !     DIVIDE ROW BY PIVOT

     KJ=K-N
     DO 170 J=1,N
        KJ=KJ+N
        !            IF (J-K) 160,170,160
        IF ( (J-K) .ne. 0 ) then
160        A(KJ) = A(KJ)/BIGA
        end if
170  end do

     !     PRODUCT OF PIVOTS

     D = D * BIGA

     !     REPLACE PIVOT BY RECIPROCAL

     !         A(KK)=1.d0/BIGA
     A(KK) = 1.0 / BIGA
180 end do

  !     FINAL ROW AND COLUMN INTERCHANGE

  K=N
190 K=(K-1)
  !      IF (K) 260,260,200
  IF ( K .gt. 0 ) then
200  I=L(K)
     !         IF (I-K) 230,230,210
     IF ( (I-K) .gt. 0 ) then
210     JQ=N*(K-1)
        JR=N*(I-1)
        DO 220 J=1,N
           JK=JQ+J
           HOLD=A(JK)
           JI=JR+J
           A(JK)=-A(JI)
           A(JI)=HOLD
220     end do
     end if
230  J=M(K)
     !         IF (J-K) 190,190,240
     IF ( (J-K) .gt. 0 ) then
240     KI=K-N
        DO 250 I=1,N
           KI=KI+N
           HOLD=A(KI)
           JI=KI-K+J
           A(KI)=-A(JI)
           A(JI)=HOLD
250     end do
     end if
     GO TO 190
  end if
260 RETURN
END SUBROUTINE MINV
!
!******************************************************************************
!
subroutine get_smacro

   USE GLOBAL
   implicit none	

      DOUBLE PRECISION :: sav6(6),sav5(5)
      DOUBLE PRECISION :: aux55(5,5),aux66(6,6),aux3333(3,3,3,3)
      DOUBLE PRECISION :: IJV(6,2)
!!
! by anand

integer ::i,j,k


ijv(1:6,1)= (/1,2,3,2,1,1/)
ijv(1:6,2)= (/1,2,3,3,3,2/)

!c
!c     OVERALL STRESS
!c
      do ii=1,3
      do jj=1,3
		  scauav(ii,jj)=0.
		  do k=1,npts3
		  do j=1,npts2
		  do i=1,npts1
			  scauav(ii,jj)=scauav(ii,jj)+sg(ii,jj,i,j,k)*wgt
		  enddo
		  enddo
		  enddo
      enddo
      enddo
!cw
!cw      enddo
!cw      enddo
!cc
!cc    MIXED BC
!cc
  do i=1,6
	  ii=ijv(i,1)
      jj=ijv(i,2)
      ddefgradmacro(ii,jj)=0.
      if(idsim(i).eq.0) then
		   do k=1,6
			   kk=ijv(k,1)
			   ll=ijv(k,2)
			   ddefgradmacro(ii,jj)=ddefgradmacro(ii,jj)+ s0(ii,jj,kk,ll)*iscau(k)*(scauchy(kk,ll)-scauav(kk,ll))
		   enddo
      endif
   enddo

   do ii=1,3
	   do jj=1,3
			ddefgradmacroacum(ii,jj)=ddefgradmacroacum(ii,jj)+ddefgradmacro(ii,jj)
	   enddo
   enddo

      write(*,*) 'DDEFGRADMACRO(1,1),(2,2)=',  ddefgradmacro(1,1),ddefgradmacro(2,2)
      call chg_basis(sav6,scauav,aux66,aux3333,2,6)

      do ii=1,5
      sav5(ii)=sav6(ii)
      enddo

      call chg_basis(sav5,sdeviat,aux55,aux3333,1,5)

      svm=0.
      do i=1,3
      do j=1,3
      svm=svm+sdeviat(i,j)*sdeviat(i,j)
      enddo
      enddo
      svm=sqrt(3./2.*svm)
 return
end
!
!*******************************************************************************************
!
      subroutine update_orient

!CC      include 'fft.dim'

   USE GLOBAL

   implicit none

      DOUBLE PRECISION :: aa(3,3),distor(3,3)
      DOUBLE PRECISION :: dnsa(3),dbsa(3)
      DOUBLE PRECISION :: rotslip(3,3),rotloc(3,3),rot(3,3)

	  ! added by anand

	  integer :: i,j,k,IPH, IS
	  DOUBLE PRECISION :: RSLBAR, RLCBAR

      RSLBAR=0.
      RLCBAR=0.

 ! MASTER DO

 do k=1,npts3
 do j=1,npts2
 do i=1,npts1

      iph=jphase(i,j,k)

  if(igas(iph).eq.0) then

!     LOCAL ROTATION RATE: ANTISYM(VELGRAD)

      do ii=1,3
      do jj=1,3
      rotloc(ii,jj)=(velgrad(ii,jj,i,j,k)-velgrad(jj,ii,i,j,k))/2.
      end do
	  end do 
 

!CCc     SLIP ROTATION RATE
!CCc
      do ii=1,3
      do jj=1,3
        aa(ii,jj)=ag(ii,jj,i,j,k)
        distor(ii,jj)=0.
      enddo
      enddo
!CCc
     do is=1,nsyst(iph)
        do ii=1,3
          dnsa(ii)=0.
          dbsa(ii)=0.
          do jj=1,3
            dnsa(ii)=dnsa(ii)+aa(ii,jj)*dnca(jj,is,iph)
            dbsa(ii)=dbsa(ii)+aa(ii,jj)*dbca(jj,is,iph)
          enddo
        enddo

        do ii=1,3
        do jj=1,3
         distor(ii,jj)=distor(ii,jj)+dbsa(ii)*dnsa(jj)*gamdot(is,i,j,k)
        enddo
        enddo
      enddo
!CCC
      do ii=1,3
      do jj=1,3
        rotslip(ii,jj)=(distor(ii,jj)-distor(jj,ii))/2.
      enddo
      enddo

!CCC     AVERAGE ROTATION RATE

      rslbar=rslbar+sqrt(rotslip(3,2)**2+rotslip(1,3)**2+rotslip(2,1)**2)*wgt
      rlcbar=rlcbar+sqrt(rotloc(3,2)**2+rotloc(1,3)**2+rotloc(2,1)**2)*wgt

!CCc    TOTAL ROTATION

     do ii=1,3
     do jj=1,3
        rot(ii,jj)=(tomtot(ii,jj)+rotloc(ii,jj)-rotslip(ii,jj))*tdot
      enddo
     enddo

!CCc     REORIENTATION

      call orient(aa,rot)

!CCc     UPDATE ORIENTATION MATRIX

     do ii=1,3
     do jj=1,3
       ag(ii,jj,i,j,k)=aa(ii,jj)
     enddo
     enddo

  endif   !  igas endif

  end do
  end do
  end do 

  write(*,*)
  write(*,*) 'AVERAGE PLASTIC ROTATION =',rslbar
  write(*,*) 'AVERAGE LOCAL ROTATION =',rlcbar
  write(*,*)

  RETURN
  END
!
!
!*******************************************************************************************
!
!*******************************************************************************************
!
      subroutine twin_Reorient

!CC      include 'fft.dim'

   USE GLOBAL

   implicit none

      DOUBLE PRECISION :: aa(3,3),distor(3,3)
      DOUBLE PRECISION :: dnsa(3),dbsa(3)
      DOUBLE PRECISION :: rotslip(3,3),rotloc(3,3),rot(3,3)

	  ! added by anand

	  integer :: i,j,k,IPH, IS
	  DOUBLE PRECISION :: TWIN_aa (3,3)

 ! MASTER DO

TWIN_aa(1,1)= -1.
TWIN_aa(1,2)=  0.
TWIN_aa(1,3)=  0.
TWIN_aa(2,1)=  0. 
TWIN_aa(2,2)=  0.064325
TWIN_aa(2,3)= -0.997929 
TWIN_aa(3,1)=  0. 
TWIN_aa(3,2)= -0.997929
TWIN_aa(3,3)=  0.064325



 do k=1,npts3
 do j=1,npts2
 do i=1,npts1

      iph=jphase(i,j,k)

if (jgrain(i,j,k) == 3) then
!CCc     UPDATE ORIENTATION MATRIX

     do ii=1,3
     do jj=1,3
       ag(ii,jj,i,j,k)=TWIN_aa(ii,jj)
     enddo
     enddo

 endif   !  igas endif

  end do
  end do
  end do 

  write(*,*)
!  write(*,*) 'AVERAGE PLASTIC ROTATION =',rslbar
!  write(*,*) 'AVERAGE LOCAL ROTATION =',rlcbar
  write(*,*)

  RETURN
  END
!
!
!***********************************************************
 subroutine orient(a,c)
 
      DOUBLE PRECISION :: a(3,3),c(3,3),th2(3,3),v(3),vbar(3)
      DOUBLE PRECISION :: th(3,3)
      DOUBLE PRECISION :: rot(3,3),anew(3,3)
      DOUBLE PRECISION :: snorm,snorm1
      INTEGER :: i,j,k

!     BUILD ROTATION TENSOR BASED ON RODRIGUES FORMULA

      v(1)=c(3,2)
      v(2)=c(1,3)
      v(3)=c(2,1)
      snorm=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
      snorm1=tan(snorm/2.)
!      if(snorm.gt.1.e-06) go to 97
!      snorm=1.
		if (snorm .le. 1.e-06)snorm = 1.


!97    do 20 i=1,3
	  do i =1,3
	      vbar(i)=snorm1*v(i)/snorm
      end do

      snorm=vbar(1)*vbar(1)+vbar(2)*vbar(2)+vbar(3)*vbar(3)
      th(3,2)=vbar(1)
      th(1,3)=vbar(2)
      th(2,1)=vbar(3)
      th(2,3)=-vbar(1)
      th(3,1)=-vbar(2)
      th(1,2)=-vbar(3)
      do i=1,3
	    th(i,i)=0.
	  end do
    
	  do i=1,3
		  do j=1,3
		  th2(i,j)=0.
			  do k=1,3
				th2(i,j)=th2(i,j)+th(i,k)*th(k,j)
			  end do
		  end do
	  end do
      
	  do i=1,3
      do j=1,3
	    rot(i,j)=(i/j)*(j/i)+2.*(th(i,j)+th2(i,j))/(1.+snorm)
      end do
	  end do

      do i=1,3
      do j=1,3
		anew(i,j)=0.
      do  k=1,3
		anew(i,j)=anew(i,j)+rot(i,k)*a(k,j)
	  end do	
      end do
	  end do

      do i=1,3
      do j=1,3
        a(i,j)=anew(i,j)
      end do
	  end do

      return
 end
!
!********************************************************************
!
!CCC *****************************************************************************
!CCC     SUBROUTINE UPDATE_CRSS_VOCE     --->      VERSION OF 23/APR/00
!CCC
!CCC     A VOCE LAW FUNCTION OF THE A!CCUMULATED SHEAR IN EACH GRAIN IS ADDED
!CCC     AS A MULTIPLICATIVE FACTOR THAT MODULATES THE ORIGINAL LINEAR HARDENING.
!CCC     THE UNITS (AND THE STRENGTH) OF THE HARDENING ARE CARRIED BY THE
!CCC     MULTIPLICATIVE FACTOR 'VOCE' (THE SAME FOR EVERY MODE).
!CCC     THE SELF & LATENT COUPLING COEFFICIENTS 'HARD' ARE DIMENSIONLESS
!CCC     CONSTANTS RELATIVE TO THE FACTOR 'VOCE'.
!CC*******************************************************************************
!CC
      SUBROUTINE HARDEN

!CC      INCLUDE 'fft.dim'

USE GLOBAL 

implicit none
! anand
integer :: i,j,k, iph, IS, JS
DOUBLE PRECISION :: GAMTOTX, DELTGAM, DTAU
!REAL :: VOCE, FACT, EXPINI, EXPDEL

      do k=1,npts3
      do j=1,npts2
      do i=1,npts1

      iph=jphase(i,j,k)

       if(igas(iph).eq.0) then

   !       GACUMGR(I,J,K)=GACUMGRaux(I,J,K)
		  GAMTOTX=GACUMGR(I,J,K)
		  
          DELTGAM=0.0
          DO IS=1,NSYST(IPH)
            DELTGAM=DELTGAM+ABS(GAMDOT(IS,I,J,K))*TDOT
          ENDDO

          DO IS=1,NSYST(IPH)
            DTAU=0.
            DO JS=1,NSYST(IPH)
              DTAU=DTAU+HARD(IS,JS,IPH)*ABS(GAMDOT(JS,I,J,K))*TDOT
            ENDDO
            TAU0 =TAU (IS,1,IPH)
            TAU1 =TAU (IS,3,IPH)  ! modified by anand from 0 index to 1
            THET0=THET(IS,1,IPH)
            THET1=THET(IS,2,IPH)
            TINY=1.E-4*TAU0

            VOCE=0.0
            IF(ABS(THET0).GT.TINY) THEN
              VOCE=THET1*DELTGAM
              IF(ABS(TAU1).GT.TINY) THEN
                FACT=ABS(THET0/TAU1)
                EXPINI=EXP(-GAMTOTX*FACT)
                EXPDEL=EXP(-DELTGAM*FACT)
                VOCE  =VOCE-(FACT*TAU1-THET1)/FACT*EXPINI*(EXPDEL-1.)-THET1/FACT*EXPINI* &
     &            (EXPDEL*((GAMTOTX+DELTGAM)*FACT+1.)-(GAMTOTX*FACT+1.))
              ENDIF
            ENDIF
            CRSS(IS,1,I,J,K)=CRSS(IS,1,I,J,K)+DTAU*VOCE/DELTGAM
            CRSS(IS,2,I,J,K)=CRSS(IS,2,I,J,K)+DTAU*VOCE/DELTGAM
			
			TRIALTAU(IS,1,I,J,K)=CRSS(IS,1,I,J,K)
			TRIALTAU(IS,2,I,J,K)=CRSS(IS,2,I,J,K)

          ENDDO
          GACUMGR(I,J,K)=GAMTOTX+DELTGAM
!CCch
!CCch        DO IS=1,NSYST(IPH)
!CCchc
!CCch          dtau=0.
!CCch            do js=1,nsyst(iph)
!CCch              dtau=dtau+hard(is,js,iph)*abs(gamdot(js,i,j,k))*tdot
!CCch            enddo
!CCch
!CCch          thet0=thet(is,0,iph)
!CCch          thet1=thet(is,1,iph)
!CCch          fact =gacumgr(i,j,k)*thet0/tau(is,3,iph)
!CCch          voce=thet1+(thet0-thet1+thet1*fact)*exp(-fact)
!CCchc
!CCch          crss(is,1,i,j,k)=crss(is,1,i,j,k)+dtau*voce
!CCch          crss(is,2,i,j,k)=crss(is,2,i,j,k)+dtau*voce
!CCchc
!CCch        ENDDO
!CCch
!CCc
      endif   !  igas endif
!CCc
      enddo
      enddo
      enddo

      END
!CCevp
!
!***************************************************************************
!
 SUBROUTINE DATA_CRYSTAL_ELAST(IPH)

      USE GLOBAL

	    double precision :: dde(3,3),xid4(3,3,3,3)
        double precision :: cc66v(6,6)
        INTEGER :: i,j,k,iph,iso
		double precision :: young, TMU,TNU,TLA
!c       UNITARY TENSORS
!c
	do i=1,3
	do j=1,3
		dde(i,j)=0.d0
		if(i.eq.j) dde(i,j)=1.d0
	enddo
	enddo

	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
	xid4(i,j,k,l)=(dde(i,k)*dde(j,l)+dde(i,l)*dde(j,k))/2.d0
	enddo
	enddo
	enddo
	enddo

	read(UR1,*) iso

	if(iso.eq.0) then
		do i=1,6
		read(UR1,*)(cc66v(i,j),j=1,6)
		enddo
		call voigt(cc66v,cc,1)
	else
		read(UR1,*) young,tnu
        tmu=young/(2.*(1.+tnu))
		tla=2.d0*tmu*tnu/(1.d0-2.d0*tnu)
		!c       bulk=2.d0*tmu*(1+tnu)/3.d0/(1.d0-2d0*tnu)
		!c
		!c	rho=3.d0*(1.d0-2.d0*tnu)/(1.d0+tnu)
		!c
	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
	cc(i,j,k,l)=tla*dde(i,j)*dde(k,l)+2.d0*tmu*xid4(i,j,k,l)
	enddo
	enddo
	enddo
	enddo

    Endif

RETURN
END

!
!*******************************************************************************
!
      subroutine evpal

      USE GLOBAL

      double precision :: xlambda(3,3),xlambda6(6)

      double precision :: sgaux(3,3),sg6(6),sg6old(6)

!cx      dimension dsg(3,3),dse(6)

      double precision :: eptaux(3,3),ept6(6)
      double precision :: edotpaux(3,3),edotp6(6)

      double precision :: dedotp66(6,6)

      double precision :: defgradaux(3,3),defgrad6(6)

      double precision :: aux66(6,6),aux3333(3,3,3,3)

      double precision :: sg66(6,6)

      double precision :: defgradceq(3,3),defgradceq6(6)
      
      double precision :: res(6)
      double precision :: xjacobinv(6,6)

! by anand

      double precision :: rss1(NSYSMX),rss2(NSYSMX)
      double precision :: sc(5,NSYSMX),rss(NSYSMX)
      double precision :: char_shear, shearstrain(3,3)
      integer :: i,j,k,itmaxal,iter1
      double precision :: sgnorm, dgnorm, erroral,erral,errald
      double precision :: dsgnorm, dsgnorm1,dsgnorm2,ddgnorm

      ERRE=0.
      ERRS=0.

   char_shear  = 0.129/n_inc   !2000.
    
	shearstrain(1,1)=0.
	shearstrain(1,2)=0.
	shearstrain(1,3)=0.
	shearstrain(2,1)=0.
	shearstrain(2,2)=-0.4989645
	shearstrain(2,3)=-0.0321624
	shearstrain(3,1)=0.
	shearstrain(3,2)=shearstrain(2,3)
	shearstrain(3,3)=0.4989645

	shearstrain = char_shear*shearstrain



   DO  k=1,npts3
   DO  j=1,npts2
   DO  i=1,npts1

   
!
! addition of the charateristic twin shear to the total strain at time t 
!   

  if (imicro .ge. 2 .and. jgrain(i,j,k) == 3  ) then
    ept(1:3,1:3,i,j,k) =  ept(1:3,1:3,i,j,k)+shearstrain(1:3,1:3)
  end if 


       do ii=1,6
       do jj=1,6
        sg66(ii,jj)=cg66(ii,jj,i,j,k)
       enddo
       enddo
      call lu_inverse(sg66,6)

       jph=jphase(i,j,k)

       do ii=1,3
       do jj=1,3
        xlambda(ii,jj)=sg(ii,jj,i,j,k)
        sgaux(ii,jj)=sg(ii,jj,i,j,k)
        eptaux(ii,jj)=ept(ii,jj,i,j,k)
        defgradaux(ii,jj)=defgrad(ii,jj,i,j,k)
       enddo
       enddo

      call chg_basis(xlambda6,xlambda,aux66,aux3333,2,6)

      call chg_basis(sg6,sgaux,aux66,aux3333,2,6)
      call chg_basis(ept6,eptaux,aux66,aux3333,2,6)
      call chg_basis(defgrad6,defgradaux,aux66,aux3333,2,6)

      sgnorm=0.
      dgnorm=0.

      do ii=1,3
      do jj=1,3
       sgnorm=sgnorm+xlambda(ii,jj)**2
       dgnorm=dgnorm+defgradaux(ii,jj)**2
      enddo
      enddo

      sgnorm=sqrt(sgnorm)
      dgnorm=sqrt(dgnorm)
!c
!cx      write(*,*) 'SGNORM,DGNORM=',sgnorm,dgnorm
!cx      pause
!c
!a      erroral=0.0000001
!a      itmaxal=100
  
      erroral=0.0000001
      itmaxal=250


      iter1=0
!cxx      erral=2*error
      erral=2*erroral
!c
!cxx      do while(iter1.lt.itmax.and.erral.gt.error)
      do while(iter1.lt.itmaxal.and.erral.gt.erroral)
      iter1=iter1+1

      sg6old=sg6

      call edotp_sg_eval(sg6,edotp6,dedotp66,i,j,k,jph)
!c
      do ii=1,6
       defgradceq6(ii)=ept6(ii)+edotp6(ii)*tdot
!cxx
!cxx    ELASTIC
!cxx
!cxx       defgradceq6(ii)=0.
!cxx
      do jj=1,6
       defgradceq6(ii)=defgradceq6(ii)+sg66(ii,jj)*sg6(jj)
      enddo
      enddo

!cx      write(*,*) 'defgrad6=',defgrad6
!cx      write(*,*) 'defgradceq6=',defgradceq6
!cx      pause

      call chg_basis(defgradceq6,defgradceq,aux66,aux3333,1,6)

!cx      if(i.eq.2.and.j.eq.2.and.k.eq.2) then
!cx      write(*,*) 'defgradceq='
!cx      do ii=1,3
!cx      write(*,*) (defgradceq(ii,jj),jj=1,3)
!cx      enddo
!cx      pause
!cx      endif

      do ii=1,6
       res(ii)=sg6(ii)-xlambda6(ii)
      do jj=1,6
       res(ii)=res(ii)+c066(ii,jj)*(defgradceq6(jj)-defgrad6(jj))
      enddo
      enddo

!cx      write(*,*) 'res=',res
!cx      pause

      do ii=1,6
      do jj=1,6
       xjacobinv(ii,jj)=(ii/jj)*(jj/ii)
      do kk=1,6
       xjacobinv(ii,jj)=xjacobinv(ii,jj)+c066(ii,kk)*(sg66(kk,jj)+dedotp66(kk,jj)*tdot)
!cxx
!cxx    ELASTIC        c066(ii,kk)*sg66(kk,jj)
      enddo
      enddo
      enddo
!c
      call lu_inverse(xjacobinv,6)
!c
      do ii=1,6
      do jj=1,6
      sg6(ii)=sg6(ii)-xjacobinv(ii,jj)*res(jj)
      enddo
      enddo
!c
      dsgnorm1=0.
      dsgnorm2=0.
      ddgnorm=0.

      do ii=1,6
       dsgnorm1=dsgnorm1+(sg6(ii)-sg6old(ii))**2
       dsgnorm2=dsgnorm2+(sg6(ii)-xlambda6(ii))**2
       ddgnorm=ddgnorm+(defgradceq6(ii)-defgrad6(ii))**2
      enddo
!c
      erral=dsgnorm1/sgnorm
      errald=ddgnorm/dgnorm

 IF(IUPHARD .EQ. 1) then


!calculate an estimate of tau_C at t+dt

!     DO is=1,nsyst(jph)
!        nrsx(is)=nrs(is,jph)
!        trialtau(is,1,i,j,k)=crss(is,1,i,j,k)
!        trialtau(is,2,i,j,k)=crss(is,2,i,j,k)
!        DO jj=1,5
!          sc(jj,is)=sch(jj,is,i,j,k)
!        ENDDO
!      ENDDO
!
!      do is=1,nsyst(jph)
!        rss(is)=sc(1,is)*sg6(1)+sc(2,is)*sg6(2)+sc(3,is)*sg6(3)+sc(4,is)*sg6(4)+sc(5,is)*sg6(5)
!        isign=1
!        if(rss(is).lt.0.) isign=2
!        rss(is)=rss(is)/trialtau(is,isign,i,j,k)
!        rss1(is)=gamd0(is,jph)*nrsx(is)*abs(rss(is)**(nrsx(is)-1))/trialtau(is,isign,i,j,k)
!        rss2(is)= gamd0(is,jph)*abs(rss(is)**(nrsx(is)))*sign(1.,rss(is))
!		GAMDOT(Is,i,j,k)=rss2(is)
!      enddo
!
!!      DO II=1,5
!!      edotp6(II)=0.
!        DO K1=1,NSYST(jph)
!          edotp6(II)=edotp6(II)+SC(II,K1)*RSS2(K1)
!        ENDDO
!      ENDDO

  
!         DELTGAM=0.0
!         DO IS=1,NSYST(JPH)
!           DELTGAM=DELTGAM+ABS(GAMDOT(IS,I,J,K))*TDOT
!         ENDDO
!
!          DO IS=1,NSYST(JPH)
!            DTAU=0.
!            DO JS=1,NSYST(JPH)
!              DTAU=DTAU+HARD(IS,JS,JPH)*ABS(GAMDOT(JS,I,J,K))*TDOT
!            ENDDO
!            TAU0 =TAU (IS,1,JPH)
!            TAU1 =TAU (IS,3,JPH)  ! modified by anand from 0 index to 1
!            THET0=THET(IS,1,JPH)
!            THET1=THET(IS,2,JPH)
!            TINY=1.E-4*TAU0
!
!            VOCE=0.0
!            IF(ABS(THET0).GT.TINY) THEN
!              VOCE=THET1*DELTGAM
!              IF(ABS(TAU1).GT.TINY) THEN
!                FACT=ABS(THET0/TAU1)
!                EXPINI=EXP(-GAMTOTX*FACT)
!                EXPDEL=EXP(-DELTGAM*FACT)
!                VOCE  =VOCE-(FACT*TAU1-THET1)/FACT*EXPINI*(EXPDEL-1.)-THET1/FACT*EXPINI* &
!     &            (EXPDEL*((GAMTOTX+DELTGAM)*FACT+1.)-(GAMTOTX*FACT+1.))
!              ENDIF
!            ENDIF
!      !      CRSS(IS,1,I,J,K)=CRSS(IS,1,I,J,K)+DTAU*VOCE/DELTGAM
!      !      CRSS(IS,2,I,J,K)=CRSS(IS,2,I,J,K)+DTAU*VOCE/DELTGAM
!			
!	 		TRIALTAU(IS,1,I,J,K)=CRSS(IS,1,I,J,K)+DTAU*VOCE/DELTGAM
!			TRIALTAU(IS,2,I,J,K)=CRSS(IS,2,I,J,K)+DTAU*VOCE/DELTGAM

!          ENDDO
     !     GACUMGR(I,J,K)=GAMTOTX+DELTGAM

	end if ! iuphard


    enddo  ! end do while (iter1)

      call chg_basis(sg6,sgaux,aux66,aux3333,1,6)
      call chg_basis(edotp6,edotpaux,aux66,aux3333,1,6)

       do ii=1,3
       do jj=1,3
        sg(ii,jj,i,j,k)=sgaux(ii,jj)  
        edotp(ii,jj,i,j,k)=edotpaux(ii,jj)
        !write(600,*)i,j,k,sg(ii,jj,i,j,k) 
       enddo
       enddo

       ERRS=ERRS+dsgnorm2*WGT
       ERRE=ERRE+ddgnorm*WGT

	END DO
	END DO
	END DO
       
      write(*,*) 'ERRE=',erre
      write(*,*) 'ERRS=',errs

      return
      end
!
!******************************************************************************
!
      subroutine edotp_sg_eval(sg6,edotp6,dedotp66,i,j,k,temjph)

    USE GLOBAL
	IMPLICIT NONE
!
      DOUBLE PRECISION ::  sg6(6)
      DOUBLE PRECISION :: edotp6(6),dedotp66(6,6)
      DOUBLE PRECISION :: rss(NSYSMX)
      DOUBLE PRECISION :: rss1(NSYSMX),rss2(NSYSMX)
      DOUBLE PRECISION :: sc(5,NSYSMX),taux(NSYSMX,2)

! added by anand
integer :: IS, JS,IPH, ISIGN,i,j,k ,temJPH, k1
!

 Jph = temJPH

      DO is=1,nsyst(jph)
        nrsx(is)=nrs(is,jph)
        taux(is,1)=TRIALTAU(IS,1,I,J,K)!crss(is,1,i,j,k)
        taux(is,2)=TRIALTAU(IS,2,I,J,K)!crss(is,2,i,j,k)
        DO jj=1,5
          sc(jj,is)=sch(jj,is,i,j,k)
        ENDDO
      ENDDO

!C     GET RESOLVED SHEAR STRESSES 'rss' AND SHEAR RATES 'gamdot'.
!C     SIGN(GAMDOT)=SIGN(RSS).
!C     NRS CAN BE EVEN OR ODD.
!C     RSS1 IS ALWAYS > 0 AND IS USED TO CALCULATE VISCOUS COMPLIANCE.

      do is=1,nsyst(jph)
        rss(is)=sc(1,is)*sg6(1)+sc(2,is)*sg6(2)+sc(3,is)*sg6(3)+sc(4,is)*sg6(4)+sc(5,is)*sg6(5)
        isign=1
        if(rss(is).lt.0.) isign=2
        rss(is)=rss(is)/taux(is,isign)
        rss1(is)=gamd0(is,jph)*nrsx(is)*abs(rss(is)**(nrsx(is)-1))/taux(is,isign)
        rss2(is)= gamd0(is,jph)*abs(rss(is)**(nrsx(is)))*sign(1.d0,rss(is))
		GAMDOT(Is,i,j,k)=rss2(is)
      enddo

      DO II=1,5
      edotp6(II)=0.
        DO K1=1,NSYST(jph)
          edotp6(II)=edotp6(II)+SC(II,K1)*RSS2(K1)
        ENDDO
      ENDDO

      edotp6(6)=0.

      DO II=1,5
      DO JJ=1,5
       dedotp66(II,JJ)=0.
       DO K1=1,NSYST(jph)
         dedotp66(II,JJ)=dedotp66(II,JJ)+SC(II,K1)*SC(JJ,K1)*RSS1(K1)
       ENDDO
      ENDDO
      ENDDO



      do ii=1,6
       dedotp66(II,6)=0.
       dedotp66(6,II)=0.
      enddo

      RETURN
   END

!************************************************
 function vm(dtensor)

	  implicit none
!c
!c     CALCULATES THE VM EQUIVALENT OF A NON-SYMMETRIC, NON-TRACELESS 
!c     (DEFORMATION GRADIENT OR VELOCITY GRADIENT) TENSOR
!c
      DOUBLE PRECISION :: dtensor(3,3),dt(3,3)
! added by anand

integer :: i,j
DOUBLE PRECISION :: trace,vm

      trace=dtensor(1,1)+dtensor(2,2)+dtensor(3,3)

      do i=1,3
      do j=1,3
       dt(i,j)=(dtensor(i,j)+dtensor(j,i))/2.-(i/j)*(j/i)*trace/3.
      enddo
      enddo

      vm=0.
      do i=1,3
      do j=1,3
      vm=vm+dt(i,j)**2
      enddo
      enddo
      vm=sqrt(2./3.*vm)

      return
end
!********************************

subroutine output_write(ioption)


use GLOBAL, only : imicro , npts1,npts2,npts3, CI ,filename,output_prefix

implicit none 

integer  :: ioption 

IF(IOPTION == 0 ) THEN 

   filename = trim(output_prefix)//'_vm.out'
   open(55,file=filename,status='unknown') 
   filename = trim(output_prefix)//'_str_str.out'
   open(56,file=filename,status='unknown')
   filename = trim(output_prefix)//'_err.out'
   open(21,file=filename,status='unknown')
   filename = trim(output_prefix)//'_conv.out'
   open(25,file=filename,status='unknown')
   filename = trim(output_prefix)//'_.Parent'
   open(510,file=filename,status ='unknown')
   filename = trim(output_prefix)//'_.Twin'
   open(511,file=filename,status ='unknown')

ELSEIF(IOPTION == 1) THEN

	write(CI,fmt='(I4)')IMICRO
	CI = ADJUSTL(CI)
	filename=TRIM(output_prefix)//'_d_'//TRIM(CI)//'.vtk' !srate
	open(91,file=filename,status='unknown') 
	filename=TRIM(output_prefix)//'_s_'//TRIM(CI)//'.vtk' !stress
	open(92,file=filename,status='unknown')
	filename=TRIM(output_prefix)//'_t_'//TRIM(CI)//'.vtk' !taufield 
	open(93,file=filename,status='unknown')
	filename=TRIM(output_prefix)//'_e_'//TRIM(CI)//'.vtk' ! strain
	open(94,file=filename,status='unknown')
	filename=TRIM(output_prefix)//'_in_'//TRIM(CI)//'.vtk' ! strain
	open(95,file=filename,status='unknown')
	filename=TRIM(output_prefix)//'_1T_'//TRIM(CI)//'.vtk' ! strain
	open(96,file=filename,status='unknown')
	filename=TRIM(output_prefix)//'_2T_'//TRIM(CI)//'.vtk' ! strain
	open(97,file=filename,status='unknown')
	filename=TRIM(output_prefix)//'_EL_'//TRIM(CI)//'.vtk' ! strain
	open(98,file=filename,status='unknown')
	filename=TRIM(output_prefix)//'_Sl_'//TRIM(CI)//'.vtk' ! strain
	open(99,file=filename,status='unknown')

!c    adding the header lines for the paraview visualization

	write(91,fmt='(A26)')'# vtk DataFile Version 2.0'
	write(91,fmt='(A2)')'# '
	write(91,fmt='(A5)')'ASCII'
	write(91,fmt='(A25)')'DATASET STRUCTURED_POINTS'
	write(91,fmt='(A11,3I5)')'DIMENSIONS ',2,npts2+1,npts3+1 !npts1+1,npts2+1,npts3+1
	write(91,fmt='(A18)')'ASPECT_RATIO 1 1 1'
	write(91,fmt='(A12)')'ORIGIN 0 0 0'
	write(91,fmt='(A10,I8)')'CELL_DATA ', 1*npts2*npts3 !npts1*npts2*npts3

	write(92,fmt='(A26)')'# vtk DataFile Version 2.0'
	write(92,fmt='(A2)')'# '
	write(92,fmt='(A5)')'ASCII'
	write(92,fmt='(A25)')'DATASET STRUCTURED_POINTS'
	write(92,fmt='(A11,3I5)')'DIMENSIONS ',2,npts2+1,npts3+1 !npts1+1,npts2+1,npts3+1
	write(92,fmt='(A18)')'ASPECT_RATIO 1 1 1'
	write(92,fmt='(A12)')'ORIGIN 0 0 0'
	write(92,fmt='(A10,I8)')'CELL_DATA ', 1*npts2*npts3 !npts1*npts2*npts3

	write(93,fmt='(A26)')'# vtk DataFile Version 2.0'
	write(93,fmt='(A2)')'# '
	write(93,fmt='(A5)')'ASCII'
	write(93,fmt='(A25)')'DATASET STRUCTURED_POINTS'
	write(93,fmt='(A11,3I5)')'DIMENSIONS ',2,npts2+1,npts3+1 !npts1+1,npts2+1,npts3+1
	write(93,fmt='(A18)')'ASPECT_RATIO 1 1 1'
	write(93,fmt='(A12)')'ORIGIN 0 0 0'
	write(93,fmt='(A10,I8)')'CELL_DATA ', 1*npts2*npts3 !npts1*npts2*npts3

	write(94,fmt='(A26)')'# vtk DataFile Version 2.0'
	write(94,fmt='(A2)')'# '
	write(94,fmt='(A5)')'ASCII'
	write(94,fmt='(A25)')'DATASET STRUCTURED_POINTS'
	write(94,fmt='(A11,3I5)')'DIMENSIONS ',2,npts2+1,npts3+1 !npts1+1,npts2+1,npts3+1
	write(94,fmt='(A18)')'ASPECT_RATIO 1 1 1'
	write(94,fmt='(A12)')'ORIGIN 0 0 0'
	write(94,fmt='(A10,I8)')'CELL_DATA ', 1*npts2*npts3 !npts1*npts2*npts3

	write(95,fmt='(A26)')'# vtk DataFile Version 2.0'
	write(95,fmt='(A2)')'# '
	write(95,fmt='(A5)')'ASCII'
	write(95,fmt='(A25)')'DATASET STRUCTURED_POINTS'
	write(95,fmt='(A11,3I5)')'DIMENSIONS ',2,npts2+1,npts3+1 !npts1+1,npts2+1,npts3+1
	write(95,fmt='(A18)')'ASPECT_RATIO 1 1 1'
	write(95,fmt='(A12)')'ORIGIN 0 0 0'
	write(95,fmt='(A10,I8)')'CELL_DATA ', 1*npts2*npts3 !npts1*npts2*npts3

	write(96,fmt='(A26)')'# vtk DataFile Version 2.0'
	write(96,fmt='(A2)')'# '
	write(96,fmt='(A5)')'ASCII'
	write(96,fmt='(A25)')'DATASET STRUCTURED_POINTS'
	write(96,fmt='(A11,3I5)')'DIMENSIONS ',2,npts2+1,npts3+1 !npts1+1,npts2+1,npts3+1
	write(96,fmt='(A18)')'ASPECT_RATIO 1 1 1'
	write(96,fmt='(A12)')'ORIGIN 0 0 0'
	write(96,fmt='(A10,I8)')'CELL_DATA ', 1*npts2*npts3 !npts1*npts2*npts3

	write(97,fmt='(A26)')'# vtk DataFile Version 2.0'
	write(97,fmt='(A2)')'# '
	write(97,fmt='(A5)')'ASCII'
	write(97,fmt='(A25)')'DATASET STRUCTURED_POINTS'
	write(97,fmt='(A11,3I5)')'DIMENSIONS ',2,npts2+1,npts3+1 !npts1+1,npts2+1,npts3+1
	write(97,fmt='(A18)')'ASPECT_RATIO 1 1 1'
	write(97,fmt='(A12)')'ORIGIN 0 0 0'
	write(97,fmt='(A10,I8)')'CELL_DATA ',1*npts2*npts3 !npts1*npts2*npts3

	write(98,fmt='(A26)')'# vtk DataFile Version 2.0'
	write(98,fmt='(A2)')'# '
	write(98,fmt='(A5)')'ASCII'
	write(98,fmt='(A25)')'DATASET STRUCTURED_POINTS'
	write(98,fmt='(A11,3I5)')'DIMENSIONS ',2,npts2+1,npts3+1 !npts1+1,npts2+1,npts3+1
	write(98,fmt='(A18)')'ASPECT_RATIO 1 1 1'
	write(98,fmt='(A12)')'ORIGIN 0 0 0'
	write(98,fmt='(A10,I8)')'CELL_DATA ', 1*npts2*npts3 !npts1*npts2*npts3

	write(99,fmt='(A26)')'# vtk DataFile Version 2.0'
	write(99,fmt='(A2)')'# '
	write(99,fmt='(A5)')'ASCII'
	write(99,fmt='(A25)')'DATASET STRUCTURED_POINTS'
	write(99,fmt='(A11,3I5)')'DIMENSIONS ',2,npts2+1,npts3+1 !npts1+1,npts2+1,npts3+1
	write(99,fmt='(A18)')'ASPECT_RATIO 1 1 1'
	write(99,fmt='(A12)')'ORIGIN 0 0 0'
	write(99,fmt='(A10,I8)')'CELL_DATA ', 1*npts2*npts3 !npts1*npts2*npts3

ENDIF ! IOPTION
    
    end subroutine


