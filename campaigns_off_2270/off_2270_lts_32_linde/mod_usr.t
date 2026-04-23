module mod_usr   ! Outer corona MHD relaxation with solar wind 
  use mod_mhd
  use mod_pfss
  implicit none
  ! some global parameters
  logical, save :: firstusrglobaldata=.true.
  real(8), allocatable :: B_init(:,:,:,:), pbc(:), rbc(:)
  integer           :: status,readwrite,blocksize
  integer*4         :: i, j, k, ilevel, nx_ss, ny_ss
  double precision  :: k_B,miu0,mass_H,usr_grav,SRadius,rhob,Tiso,rhob1,qs
  
  ! parameters of Parker's solar wind
  double precision  :: rc, Vs
  double precision  :: Vout, V_surface
  integer           :: ix^D,nth,nph,nvec,nrr
  double precision :: dth,dph,dr0,mag_divisor
  integer             :: ixg2,ixg3
  double precision    :: xlen2, xlen3
  integer, dimension(4) :: shp, shp_bound

  ! global parameters for RBSL
  integer :: np,nnp
  double precision, allocatable :: x_axis(:,:)
  double precision :: F_flx, a0
  logical          :: I_Helix
  ! global parameters for RBSL boundary
  integer          :: ixgI^L,ixgO^L 
  integer          :: m,n,p
  double precision, allocatable :: brtp_rbsl(:,:,:,:)
  double precision, allocatable :: brtp_rbsl_inner(:,:,:,:)
  logical, save :: bound_rbsl_top_ready=.false.
  double precision :: dx01, dx02, dx03

  ! global parameters for RBSL inner boundary
  double precision, allocatable :: Bfr_bound(:,:,:,:) 
  logical, save :: bound_rbsl_ready = .false. 
  double precision :: dra1, dth1, dph1

  ! global parameters for prominence density along RBSL axis
  double precision, allocatable :: s_axis(:), t_axis(:,:), n_axis(:,:), b_axis(:,:)
  logical :: prominence_ready=.false.
  double precision :: mm2n
  double precision :: prom_Crho, prom_zc, prom_zbase, prom_Fmax
  double precision :: prom_l_s, prom_w_s, prom_l_n, prom_w_n, prom_l_b, prom_w_b

contains

  !==============================================================================
  ! Purpose: to include global parameters, set user methods, set coordinate 
  !          system and activate physics module.
  !==============================================================================
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    !usr_transform_w     => transform_w_prominence

    call set_coordinate_system("spherical")
    call mhd_activate()
  end subroutine usr_init

  !==============================================================================
  ! Purpose: to initialize user public parameters and reset global parameters.
  !          Input data are also read here.
  !==============================================================================
  subroutine initglobaldata_usr()
    use mod_global_parameters
    character(len=70) :: file_path, file_Bfr_bound, file_brtp_rbsl
    logical :: flip
    double precision, allocatable :: b_r0(:,:),theta(:),phi(:)
    integer :: xm,ym,amode,file_handle
    integer :: lat0,lat1,lon0,lon1
    integer :: iphi, itheta, unit
    integer, dimension(MPI_STATUS_SIZE) :: statuss
    character(len=100)    :: filename
    real(8)               :: dr
    integer               :: ibc

    ! normalization unit in CGS Unit
    k_B = 1.3806d-16          ! erg*K^-1
    miu0 = 4.d0*dpi           ! Gauss^2 cm^2 dyne^-1
    mass_H = 1.67262d-24      ! g
    unit_length        = 6.955d10 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3
    unit_density       = 1.4d0*mass_H*unit_numberdensity               ! 2.341668000000000E-015 g*cm^-3
    unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature ! 0.317538000000000 erg*cm^-3
    unit_magneticfield = dsqrt(miu0*unit_pressure)                     ! 1.99757357615242 Gauss
    unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density)   ! 1.16448846777562E007 cm/s = 116.45 km/s
    unit_time          = unit_length/unit_velocity                     ! 5972.5794 s = 99.543 min
    
    !R_s =3.d0
    usr_grav=-2.74d4*unit_length/unit_velocity**2  ! solar gravity
    SRadius=6.955d10/unit_length                   ! Solar radius
    mag_divisor = 1.0d0
    rc = 3.45d0
    Vs = 117.54d0 
    !rhob = 5.0d9/100.0d0/unit_numberdensity 
    rhob = 1.0d8/unit_numberdensity
    Tiso = 2.01d0
    qs = dble(qstretch_baselevel(1))
    print*, 'qstretch level:', qs

    filename = './initial/OFF_combined_lmax10_q1.008_nr600.bin'
    call read_initial_magnetic_field(filename, shp, B_init)

    nrr = domain_nx1
    nth = domain_nx2
    nph = domain_nx3
    dth = dble(xprobmax2-xprobmin2)/dble(nth-1)
    dph = dble(xprobmax3-xprobmin3)/dble(nph-1)
    !dr0 = dble(xprobmax1-xprobmin1)/dble(nrr-0)
    dr0 = dble(xprobmax1-xprobmin1)*dble(1.d0-qs)/dble(1.d0-qs**dble(nrr-0))

    call cal_parker_solar_wind(1.0d0,rc,vs,V_surface)

    allocate(pbc(nghostcells))
    allocate(rbc(nghostcells))
    dr = 0
    do ibc=nghostcells,1,-1
        dr = dr+dr0*qs**(-ibc)
        call cal_parker_solar_wind(1-dr,rc,vs,Vout)
        rbc(ibc) = (rhob*V_surface)/(Vout*(1-dr)**2)
        pbc(ibc) = rbc(ibc)*Tiso
        if (mype==0) then
        print*, 'PSW at ', ibc, 'ghost layer'
        print*, 'Velocity at current layer', Vout
        print*, 'number density at current layer', rbc(ibc)
        print*, 'p_ at current layer', pbc(ibc)
        endif
    end do

    if (mype == 0) then
       call cal_parker_solar_wind(xprobmax1, rc, vs, vout)
       print *,'number density at maxmal r', (rhob*V_surface)/(Vout*xprobmax1**2) 
       print *,'velocity ar maxmal r', vout
    endif

  end subroutine initglobaldata_usr

  !==============================================================================
  ! Purpose: to read initial magnetic field
  !==============================================================================  
  subroutine read_initial_magnetic_field(filename, shp, B_init)
      use mod_global_parameters
      character(len=*), intent(in)               :: filename
      double precision, allocatable, intent(out) :: B_init(:,:,:,:)
      integer, dimension(4), intent(out)         :: shp
      integer                                    :: ios
      integer                                    :: unit
      
      open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
      read(unit) shp

      allocate(B_init(shp(1),shp(2),shp(3),shp(4)))
      read(unit) B_init
      close(unit)
      B_init = B_init(:,:,shp(3):1:-1,:)

      if(mype==0) then
      print *,'MHD simulation initialized by amrvac PFSS model'
      print *, '  Initial Data Shape: ', shp
      print *, '  Load Initial Magnetic field data from: ', filename
      print *, 'First element before normalization: ', B_init(1,1,1,1)
      print *, 'Last  element before normalization: ', B_init(size(B_init,1),size(B_init,2),size(B_init,3),size(B_init,4))
      endif

      B_init = B_init/unit_magneticfield/mag_divisor
      if(mype==0) then
      print*,'    min, max values of initial B (normalized):', minval(B_init), maxval(B_init)
      endif
  end subroutine read_initial_magnetic_field

  !==============================================================================
  ! Purpose: to initialize the initial condition
  !==============================================================================
  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    use mod_global_parameters    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: qs
    integer          :: irr,ith,iph
    logical, save:: first=.true.
    

    if (first .and. mype==0) then
      print *,'Relax a solor wind model with politropic MHD'
      print *,'User Tiso(MK): ', Tiso
      print *,'V_surface(unit_velocity): ', V_surface
      first=.false.
    end if
    
    qs = dble(qstretch_baselevel(1))
    w(ixO^S,mom(1:3))=0.0d0
    !w(ixO^S,mag(1):mag(ndir))=B_init(ixO^S,1:ndir)
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      call cal_parker_solar_wind(x(ix^D,1), rc, vs, vout)
      w(ix^D,rho_)  = (rhob*V_surface)/(Vout*x(ix^D,1)**2)
      w(ix^D,mom(1))= Vout*w(ix^D,rho_)
      irr = ceiling(log(1-(1-qs)*(x(ix^D,1)-xprobmin1)/dr0)/log(qs))
      !irr = ceiling((log(1-(1-qs)*(x(ix^D,1)-xprobmin1)/dr0)-log(1+1/(2*qs)*(1-qs)))/log(qs))
      !irr = ceiling((x(ix1,ix2,ix3,1)-xprobmin1)/dr0)
      ith = ceiling((x(ix^D,2)-xprobmin2)/dth) 
      iph = ceiling((x(ix^D,3)-xprobmin3)/dph)
      w(ix^D,mag(1:3))=B_init(1:3,irr,ith,iph)
    {end do\}
    !w(ixO^S,rho_) = rhob1*dexp(usr_grav*(SRadius**2)/Tiso*(1.d0/SRadius-1.d0/x(ixO^S,1))) ! isotermal atomosphere
    w(ixO^S,p_)   = w(ixO^S,rho_)*Tiso/(mhd_gamma-1.0d0)
    
    if(mhd_glm) w(ixO^S,psi_)=0.d0
    !call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  !==============================================================================
  ! Purpose: Calculate the Parker solar wind solution with Lambert function
  !==============================================================================
  subroutine cal_parker_solar_wind(ra_pksw, rc_pksw, vs_pksw, v_pksw)
  use mod_global_parameters
  real*8, intent(in)  :: ra_pksw, rc_pksw, vs_pksw
  real*8, intent(out) :: v_pksw
  real*8              :: W0_pksw, Wn1_pksw, Dr_pksw, nDr_pksw, plw, L_pksw, M_pksw


   Dr_pksw = ((ra_pksw/rc_pksw)**(-4.0d0))*dexp(4.0d0*(1.0d0-(rc_pksw/ra_pksw))-1.0d0)
   nDr_pksw = -1.0d0*Dr_pksw
   if (ra_pksw .le. rc_pksw-0.05d0) then
       if (nDr_pksw .gt. -0.25) then
            W0_pksw = nDr_pksw-nDr_pksw**2.0d0 + 1.5d0*nDr_pksw*nDr_pksw*nDr_pksw
       endif
       if (nDr_pksw .le. -0.25) then
           plw=sqrt(2.0d0*(exp(1.0d0)*nDr_pksw+1.0d0))
           W0_pksw=-1.0d0+plw-(1.0d0/3.0d0)*plw**2.0d0+(11.0d0/72.0d0)*plw*plw*plw
       endif
     v_pksw = dsqrt(-1.0d0*vs_pksw*vs_pksw*W0_pksw)
   endif

    if (ra_pksw .ge. rc_pksw-0.05d0 .and. ra_pksw .le. rc_pksw+0.05d0) then
     v_pksw = vs_pksw*sqrt(3.0d0-2.0d0*rc_pksw/ra_pksw)
    endif


   if (ra_pksw .gt. rc_pksw+0.05d0) then
      L_pksw=log(abs(nDr_pksw))
      M_pksw=log(abs(L_pksw))
      Wn1_pksw = L_pksw - M_pksw + (M_pksw/L_pksw) + &
                 M_pksw*(M_pksw-2.0d0)/(2.0d0*L_pksw*L_pksw)+&
                 M_pksw*(6.0d0-9.0d0*M_pksw+2.0d0*M_pksw*M_pksw)/(6.0d0*L_pksw*L_pksw*L_pksw)
      if (abs(nDr_pksw) .ge. 0.28796090d0) then
         plw=-sqrt(2.0d0*(exp(1.0d0)*nDr_pksw+1.0d0))
         Wn1_pksw=-1.0d0+plw-(1.0d0/3.0d0)*plw*plw+(11.0d0/72.0d0)*plw*plw*plw
      endif
     v_pksw = dsqrt(-1.0d0*vs_pksw*vs_pksw*Wn1_pksw)
   endif
   ! normalize the unit
     v_pksw = v_pksw*1.0d5/unit_velocity  ! Km/s -> normalized velocity

  end subroutine

  !==============================================================================
  ! Purpose: convert vectors in Cartesian coordinates to spherical ones
  !==============================================================================
  subroutine Cart2SphereVector(ixI^L,ixO^L,x,A_in,A_out)

    integer,intent(in)           :: ixI^L,ixO^L
    double precision,intent(in)  :: x(ixI^S,1:ndim)
    double precision,intent(in)  :: A_in(ixI^S,1:ndim)
    double precision,intent(out) :: A_out(ixI^S,1:ndim)

    double precision :: raddeg = dpi/ 180.
    double precision :: lon(ixI^S),lat(ixI^S)
    double precision :: bxCart(ixI^S),byCart(ixI^S),bzCart(ixI^S)
    double precision :: br(ixI^S),bth(ixI^S),bph(ixI^S)
    double precision :: a11(ixI^S),a12(ixI^S),a13(ixI^S)
    double precision :: a21(ixI^S),a22(ixI^S),a23(ixI^S)
    double precision :: a31(ixI^S),a32(ixI^S),a33(ixI^S)
    double precision :: latc(ixI^S),lonc(ixI^S),pAng(ixI^S)

    bxCart=zero
    byCart=zero
    bzCart=zero
       lon=zero
       lat=zero
       bph=zero
       bth=zero
        br=zero
     A_out=zero

    bxCart(ixO^S) = A_in(ixO^S,2)
    byCart(ixO^S) = A_in(ixO^S,3)
    bzCart(ixO^S) = A_in(ixO^S,1)
    lon(ixO^S) =  x(ixO^S,3)
    lat(ixO^S) =  0.5d0*dpi - x(ixO^S,2)
    latc = 0.0d0
    lonc = 0.0d0
    pAng = 0.0d0

    a11 = -sin(latc) * sin(pAng) * sin(lon - lonc) + cos(pAng) * cos(lon - lonc)
    a12 =  sin(latc) * cos(pAng) * sin(lon - lonc) + sin(pAng) * cos(lon - lonc)
    a13 = -cos(latc) * sin(lon - lonc)
    a21 = -sin(lat) * (sin(latc) * sin(pAng) * cos(lon - lonc) + cos(pAng) * sin(lon - lonc)) - cos(lat) * cos(latc) * sin(pAng)
    a22 =  sin(lat) * (sin(latc) * cos(pAng) * cos(lon - lonc) - sin(pAng) * sin(lon - lonc)) + cos(lat) * cos(latc) * cos(pAng)
    a23 = -cos(latc) * sin(lat) * cos(lon - lonc) + sin(latc) * cos(lat)
    a31 =  cos(lat) * (sin(latc) * sin(pAng) * cos(lon - lonc) + cos(pAng) * sin(lon - lonc)) - sin(lat) * cos(latc) * sin(pAng)
    a32 = -cos(lat) * (sin(latc) * cos(pAng) * cos(lon - lonc) - sin(pAng) * sin(lon - lonc)) + sin(lat) * cos(latc) * cos(pAng)
    a33 =  cos(lat) * cos(latc) * cos(lon - lonc) + sin(lat) * sin(latc)

    bph(ixO^S) = a11(ixO^S) * bxCart(ixO^S) +& 
                 a12(ixO^S) * byCart(ixO^S) +& 
                 a13(ixO^S) * bzCart(ixO^S)
    bth(ixO^S) = a21(ixO^S) * bxCart(ixO^S) +&
                 a22(ixO^S) * byCart(ixO^S) +&
                 a23(ixO^S) * bzCart(ixO^S)
     br(ixO^S) = a31(ixO^S) * bxCart(ixO^S) +&
                 a32(ixO^S) * byCart(ixO^S) +&
                 a33(ixO^S) * bzCart(ixO^S)

    A_out(ixO^S,1) =         br(ixO^S)
    A_out(ixO^S,2) = -1.0d0*bth(ixO^S)
    A_out(ixO^S,3) =        bph(ixO^S)

    !print*,'Finish Cart2SphereVector!'   

  end subroutine Cart2SphereVector

subroutine Spherical_curlvector_usr(qvec,ixI^L,ixO^L,xS,curlvec)
  use mod_global_parameters
  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in)  :: qvec(ixI^S,1:3), xS(ixI^S,1:3)
  double precision, intent(out) :: curlvec(ixO^S,1:3)
  integer :: ixA^L, hxO^L, jxO^L, idir, jdir, kdir
  double precision :: tmp(ixI^S), tmp2(ixI^S), epssin

  ixA^L=ixO^L^LADD1;
  if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
    call mpistop("Error in curlvector: Non-conforming input limits")
  curlvec(ixO^S,1:3)=zero

  epssin = 1.0d-6
  do idir=1,3; do jdir=1,ndim; do kdir=1,3
    if(lvc(idir,jdir,kdir)/=0)then
      tmp(ixA^S)=qvec(ixA^S,kdir)
      hxO^L=ixO^L-kr(jdir,^D);
      jxO^L=ixO^L+kr(jdir,^D);
      select case(jdir)
      case(1)
      tmp(ixA^S)=tmp(ixA^S)*xS(ixA^S,1)
      tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((xS(jxO^S,1)-xS(hxO^S,1))*xS(ixO^S,1))
      {^NOONED    case(2)
      if(idir==1) tmp(ixA^S)=tmp(ixA^S)*dsin(xS(ixA^S,2))
      tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((xS(jxO^S,2)-xS(hxO^S,2))*xS(ixO^S,1))
      if(idir==1) tmp2(ixO^S)=tmp2(ixO^S)/max(dabs(dsin(xS(ixO^S,2))),epssin)
      }
      {^IFTHREED  case(3)
      tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((xS(jxO^S,3)-xS(hxO^S,3))*xS(ixO^S,1)*dsin(xS(ixO^S,2)))
      }
      end select
      if(lvc(idir,jdir,kdir)==1)then
        curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
      else
        curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
      endif
    endif
  enddo; enddo; enddo;
end subroutine Spherical_curlvector_usr

  !==============================================================================
  ! Purpose: to provide special boundary conditions set by users.
  !==============================================================================
  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    
    double precision :: Qp(ixI^S),xlen^D
    double precision :: r_ghost, r_mirror
    double precision, allocatable :: magCC(:,:,:),xex(:,:,:)
    integer :: ix^D,ixOs^L,jxO^L,idir,ixbc^D,ixA^L
    
    !! MHD extra parameters
    integer :: ixIM^L, ith, iph, ira
    integer :: ith1, iph1
    integer :: ira2, ith2, iph2
    double precision :: tmp1(ixI^S),tmp2(ixI^S),pth(ixI^S),tmpB(ixO^S,1:ndim)
    double precision :: coeffrho
    double precision :: q,q1,a0,a1,a2,a3,b0,b1,b2,c1,c2,c3,c4
    double precision :: BfrS_interp(3)
    
    q  = 1.008d0
    q1 = 1.0d0/q
    a0 = -1.0d0-1.0d0/(1.0d0+q)-1/(1.0d0+q+q**2)
    a1 = 1.0d0+1.0d0/q**2+1.0d0/q**2
    a2 = -(1.0d0+q+q**2)/(q**3+q**4)
    a3 = 1.0d0/(q**3+q**4+q**5)
    b0 = -(2.0d0+q1)/(1.0d0+q1)
    b1 = 1.0d0+1.0d0/q1
    b2 = -1.0d0/(1.0d0+q1**2)
    ! constant value extrapolation in r-maximal boundary 
    c1 = -24.0d0*q1**6/(1.0d0-2.0d0*q1-3.0d0*q1**2+2.0d0*q1**3+8.0d0*q1**4+12.0d0*q1**5-24.0d0*q1**6)
    c2 = 6.0d0*q1**3/(2.0d0-7.0d0*q1+2.0d0*q1**2+14.0d0*q1**3-12.0d0*q1**4)
    c3 = 8.0d0*q1/(-6.0d0+17.0d0*q1+6.0d0*q1**2-51.0d0*q1**3+36.0d0*q1**4)
    c4 = 3.0d0/(3.0d0-4.0d0*q1-6.0d0*q1**2-4.0d0*q1**3+16.0d0*q1**4+24.0d0*q1**5-32.0d0*q1**6)

    nrr = domain_nx1
    nth = domain_nx2
    nph = domain_nx3
    dth = dble(xprobmax2-xprobmin2)/dble(nth-1)
    dph = dble(xprobmax3-xprobmin3)/dble(nph-1)

    !if(mhd_glm) w(ixO^S,psi_)=0.d0
    
    select case(iB)
     case(1)
       !if(mhd_glm) w(ixO^S,psi_)=0.d0
       !call mhd_to_primitive(ixI^L,IxO^L,w,x)
       w(ixO^S,mom(:))=0.0d0
       ! velocity and density are fixed to PKSW solution
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
         call cal_parker_solar_wind(x(ix^D,1), rc, vs, vout)
         w(ix^D,rho_)  = (rhob*V_surface)/(Vout*x(ix^D,1)**2)
         w(ix^D,p_)    = w(ix^D,rho_)*Tiso/(mhd_gamma-1.0d0)
       {end do\}
       w(ixO^S,mom(1)) = zero
       w(ixO^S,mom(2)) = zero
       w(ixO^S,mom(3)) = zero
       !if(mype==0) then
       !  print *, 'maxval, minval of number density:', maxval(w(ixO^S, rho_)), minval(w(ixO^S, rho_))
       !  print *, 'maxval, minval of p:', maxval(w(ixO^S, p_)), minval(w(ixO^S, p_))
       !endif
       !w(ixO^S,mom(2)) = -w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))
       !w(ixO^S,mom(3)) = -w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))
!       w(ixO^S,mom(2)) = -w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))/&
!                          w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) 
!       w(ixO^S,mom(3)) = -w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))/&
!                          w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
!       ! pressure use the equal gradient extrapolation of 3rd
!       do ix1=ixOmax1,ixOmin1,-1
!         w(ix1^%1ixO^S,p_)=(a0-a1)/a0*w(ix1+1^%1ixO^S,p_)+&
!                           (a1-a2)/a0*w(ix1+2^%1ixO^S,p_)+&
!                           (a2-a3)/a0*w(ix1+3^%1ixO^S,p_)+&
!                           (a3- 0)/a0*w(ix1+4^%1ixO^S,p_)
!       end do
!    ! pressure use the isothermal atmosphere
!       pth(ixOmax1+1^%1ixO^S) = w(ixOmax1+1^%1ixO^S,p_)/w(ixOmax1+1^%1ixO^S,rho_) ![T]
!       do ix1=ixOmax1, ixOmin1, -1
!         w(ix1^%1ixO^S,p_)   = pth(ixOmax1+1^%1ixO^S)*w(ix1^%1ixO^S,rho_)
!       end do
       ! magnetic field: first ghost cell with line-tied condition; the other use the equal gradient extrapolation
       ! w(ixOmax1^%1ixO^S,mag(1))   = w(ixOmax1+1^%1ixO^S,mag(1))
       !if (.not. bound_rbsl_ready) call bound_Brtp_rbsl()
       do ix2=ixOmin2,ixOmax2
         do ix3 = ixOmin3,ixOmax3
           ! fixed to OFF field
           ith1 = max(1, min(nth, ceiling((x(ixOmax1,ix2,ix3,2)-xprobmin2)/dth)))
           iph1 = max(1, min(nph, ceiling((x(ixOmax1,ix2,ix3,3)-xprobmin3)/dph)))
           w(ixOmax1,ix2,ix3,mag(1:3)) = B_init(1:3,1,ith1,iph1)
         end do
       end do
       !w(ixOmax1^%1ixO^S,mag(2:3)) = (a0-a1)/a0*w(ixOmax1+1^%1ixO^S,mag(2:3))+&
       !                              (a1-a2)/a0*w(ixOmax1+2^%1ixO^S,mag(2:3))+&
       !                              (a2-a3)/a0*w(ixOmax1+3^%1ixO^S,mag(2:3))+&
       !                              (a3- 0)/a0*w(ixOmax1+4^%1ixO^S,mag(2:3))
       do ix1=ixOmax1-1,ixOmin1,-1
         w(ix1^%1ixO^S,mag(1:3))=(a0-a1)/a0*w(ix1+1^%1ixO^S,mag(1:3))+&
                                 (a1-a2)/a0*w(ix1+2^%1ixO^S,mag(1:3))+&
                                 (a2-a3)/a0*w(ix1+3^%1ixO^S,mag(1:3))+&
                                 (a3- 0)/a0*w(ix1+4^%1ixO^S,mag(1:3))
       end do
       !call mhd_to_conserved(ixI^L,ixO^L,w,x)

     case(2)
       !if(mhd_glm) w(ixO^S,psi_)=0.d0
       !call mhd_to_primitive(ixI^L,ixA^L,w,x)
       do ix1 = ixOmin1,ixOmax1
         ! constant-value extrapolation
         !w(ix1^%1ixO^S,rho_) = c1*w(ix1-1^%1ixO^S,rho_)+&
         !                      c2*w(ix1-2^%1ixO^S,rho_)+&
         !                      c3*w(ix1-3^%1ixO^S,rho_)+&
         !                      c4*w(ix1-4^%1ixO^S,rho_)
         !w(ix1^%1ixO^S,p_)   = c1*w(ix1-1^%1ixO^S,p_)+&
         !                      c2*w(ix1-2^%1ixO^S,p_)+&
         !                      c3*w(ix1-3^%1ixO^S,p_)+&
         !                      c4*w(ix1-4^%1ixO^S,p_)
         !w(ix1^%1ixO^S,mom(1)) = c1*w(ix1-1^%1ixO^S,mom(1))+&
         !                        c2*w(ix1-2^%1ixO^S,mom(1))+&
         !                        c3*w(ix1-3^%1ixO^S,mom(1))+&
         !                        c4*w(ix1-4^%1ixO^S,mom(1))

         ! zero-gradient extrapolation
         w(ix1^%1ixO^S,rho_) = (-b1/b0)*w(ix1-1^%1ixO^S,rho_)+&
                               (-b2/b0)*w(ix1-2^%1ixO^S,rho_)
         w(ix1^%1ixO^S,p_)   = (-b1/b0)*w(ix1-1^%1ixO^S,p_)+&
                               (-b2/b0)*w(ix1-2^%1ixO^S,p_)
         w(ix1^%1ixO^S,mom(1)) = (-b1/b0)*w(ix1-1^%1ixO^S,mom(1))+&
                                 (-b2/b0)*w(ix1-2^%1ixO^S,mom(1))
         w(ix1^%1ixO^S,mom(2:3)) = zero
         w(ix1^%1ixO^S,mag(2:3)) = zero
         ! zero-gradient extra, r^2 Br conservation
         w(ix1^%1ixO^S,mag(1)) = ((-b1/b0)*w(ix1-1^%1ixO^S,mag(1))*x(ix1-1^%1ixO^S,1)**2+&
                                  (-b2/b0)*w(ix1-2^%1ixO^S,mag(1))*x(ix1-2^%1ixO^S,1)**2)/x(ix1^%1ixO^S,1)**2 
                                  !magnetic flux constant 
!         w(ix1^%1ixO^S,mag(:)) = c1*w(ix1-1^%1ixO^S,mag(:))+&
!                                 c2*w(ix1-2^%1ixO^S,mag(:))+&
!                                 c3*w(ix1-3^%1ixO^S,mag(:))+&
!                                 c4*w(ix1-4^%1ixO^S,mag(:))
       end do
       
       ! rho, mr fixed to PKSW solution
       !{do ix^DB=ixOmin^DB,ixOmax^DB\}
         !call cal_parker_solar_wind(x(ix^D,1), rc, vs, vout)
         !w(ix^D,rho_)  = (rhob*V_surface)/(Vout*x(ix^D,1)**2)
         !w(ix^D,mom(1))= vout*w(ix^D,rho_)
       !{end do\}
       
       ! mr no-inflow
       !do ix1 = ixOmin1,ixOmax1
       !  w(ix1^%1ixO^S,mom(1)) = w(ixOmin1-1^%1ixO^S,mom(1))
       !  where(w(ix1^%1ixO^S,mom(1))<0.d0)
       !    w(ix1^%1ixO^S,mom(1))=0.d0
       !  end where
       !end do

       ! mt, mp =0
       !w(ix1^%1ixO^S,mom(2:3)) =  0

       ! mt, mp inverse-symmetric
       !w(ixO^S,mom(2)) = -w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))
       !w(ixO^S,mom(3)) = -w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))
       
       !do ix1 = ixOmin1,ixOmax1
          ! r^2*rho is continuous
          ! w(ix1^%1ixO^S,rho_) = w(ixOmin1-1^%1ixO^S,rho_)*(x(ixOmin1-1^%1ixO^S,1)/x(ix1^%1ixO^S,1))**2
          ! v_r no-inflow is continuous
          !waii w(ix1^%1ixO^S,mom(1)) = w(ixOmin1-1^%1ixO^S,mom(1))
          !where(w(ix1^%1ixO^S,mom(1))<0.d0)
          !  w(ix1^%1ixO^S,mom(1))=0.d0
          !end where
          ! r*v_theta is continuous, angular momentum conservation
          !w(ix1^%1ixO^S,mom(2)) = w(ixOmin1-1^%1ixO^S,mom(2))*x(ixOmin1-1^%1ixO^S,1)/x(ix1^%1ixO^S,1)
          ! r*v_phi is continuous, angular momentum conservation
          !w(ix1^%1ixO^S,mom(3)) = w(ixOmin1-1^%1ixO^S,mom(3))*x(ixOmin1-1^%1ixO^S,1)/x(ix1^%1ixO^S,1)
       !end do

       ! p zero-gradient
       !do ix1=ixOmin1,ixOmax2
       !  w(ix1^%1ixO^S,p_)=(-b1/b0)*w(ix1-1^%1ixO^S,p_)+&
       !                    (-b2/b0)*w(ix1-2^%1ixO^S,p_)
       !end do

       ! p: isothermal
       !pth(ixOmin1-1^%1ixO^S) = w(ixOmin1-1^%1ixO^S,p_)/w(ixOmin1-1^%1ixO^S,rho_) ![T]
       !do ix1=ixOmin1,ixOmax1
         ! isothermal boundary condition
       !  w(ix1^%1ixO^S,p_)   = pth(ixOmin1-1^%1ixO^S)*w(ix1^%1ixO^S,rho_)
       !enddo
       
       ! magnetic field: Br use zero-gradient extrapolation; Btp fix to zero
       !do ix1=ixOmin1,ixOmax1
       !  w(ix1^%1ixO^S, mag(1))=(-b1/b0)*w(ix1-1^%1ixO^S,mag(1))+&
       !                         (-b2/b0)*w(ix1-2^%1ixO^S,mag(1))
       ! enddo
       !w(ixO^S,mag(2:3))=0.d0
       !call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case(3)

     case(4)

     case(5)

     case(6)

     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  !==============================================================================
  ! Purpose: get gravity field
  !==============================================================================
  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*SRadius**2/(x(ixO^S,1))**2
  end subroutine

  !==============================================================================
  ! Purpose: get gravity field
  !==============================================================================
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)
    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,1)=ggrid(ixO^S)
  end subroutine gravity

  !==============================================================================
  ! Purpose: Enforce additional refinement or coarsening. One can use the
  !          coordinate info in x and/or time qt=t_n and w(t_n) values w.
  !==============================================================================
    subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision :: joverb_cr, dxmin
    integer, intent(inout) :: refine, coarsen
    double precision  :: current(ixI^S,1:ndir), btotal(ixI^S,1:ndir), absj(ixO^S), absb(ixO^S)
    integer :: idirmin

    dxmin = minval(dxlevel(1:ndim))

    joverb_cr = 1.d0 / dxmin
    btotal = w(ixI^S,mag(:))
    call curlvector(btotal,ixI^L,ixO^L,current,idirmin,1,ndir)
    absj(ixO^S) = dsqrt(sum(current(ixO^S,:)**2,dim=ndim+1))
    absb(ixO^S) = dsqrt(sum(btotal(ixO^S,:)**2,dim=ndim+1))

    !write(*,'(A,I6,A,I2,A,ES12.5,A,ES12.5)') 'igrid=', igrid, ' level=', level, &
    ! ' dx_min_current=', dxmin, ' joverb_cr=', joverb_cr

     if ((level<refine_max_level) .and. any( (absj/absb)> joverb_cr )) then
        refine=1
        coarsen=-1
     else if ((level<(refine_max_level-1)) .and. any( (absj/absb)>(joverb_cr-0.3/dxmin) )) then
        refine=1
        coarsen=-1
      else if ((level<(refine_max_level-2)) .and. any( (absj/absb)>(joverb_cr-0.5/dxmin) )) then
        refine=1
        coarsen=-1
      else
        refine=-1
        coarsen=0
      end if
    
  end subroutine special_refine_grid

  !==============================================================================
  ! Purpose: 
  !   this subroutine can be used in convert, to add auxiliary variables to the
  !   converted output file, for further analysis using tecplot, paraview, ....
  !   these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  !   the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  !   corresponding normalization values (default value 1)
  !==============================================================================
  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: w1(ixI^S,nw+nwauxio)
    double precision                   :: qvec(ixI^S,1:ndim)
    double precision                   :: tmp(ixI^S)
    double precision                   :: current(ixI^S,1:ndim), Lorentz(ixI^S,1:ndim)
    double precision                   :: xlen^D
    integer :: ix^D,ixbc^D,idirmin,idir,jdir,kdir
    double precision                   :: divb(ixI^S),B2(ixI^S),dip(ixI^S)
    double precision                   :: divv(ixI^S), div_free_fi(ixI^S)
    double precision                   :: vels(ixI^S,1:ndim)
    ! parameters for magnetic field dips
    double precision :: dBr_dr(ixI^S), dBr_dth(ixI^S), dBr_dph(ixI^S)
    double precision :: dip_dir(ixI^S), r(ixI^S), th(ixI^S), sinth(ixI^S)
    double precision :: Bmag(ixI^S), epsBr(ixI^S), eps_sin

    ! output Brtp
    w(ixO^S,nw+1:nw+nwauxio)=zero
    w(ixO^S,nw+1 )=w(ixO^S,mag(1))*unit_magneticfield
    w(ixO^S,nw+2 )=w(ixO^S,mag(2))*unit_magneticfield
    w(ixO^S,nw+3 )=w(ixO^S,mag(3))*unit_magneticfield
    ! output current Jrtp
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    B2(ixO^S)=sum(w(ixO^S,mag(1:ndir))**2,dim=ndim+1)
    w(ixO^S,nw+4 )=current(ixO^S,1)
    w(ixO^S,nw+5 )=current(ixO^S,2)
    w(ixO^S,nw+6 )=current(ixO^S,3)
    ! output divb
    call get_divb(w,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+7 )=tmp(ixO^S)
    ! output thermal pressure
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+8 )=tmp(ixO^S)
    ! output radial velocity Vr
    w(ixO^S,nw+9 )=w(ixO^S,mom(1))/w(ixO^S,rho_)*unit_velocity/1.d5
    vels(ixO^S, 1)=w(ixO^S,mom(1))/w(ixO^S,rho_)
    vels(ixO^S, 2)=w(ixO^S,mom(2))/w(ixO^S,rho_)
    vels(ixO^S, 3)=w(ixO^S,mom(3))/w(ixO^S,rho_)
    ! output divv
    !call divvector(vels(ixI^S,1:ndir),ixI^L,ixO^L,divv)
    !w(ixO^S,nw+10)=divv(ixO^S)
    !call get_normalized_divb(w,ixI^L,ixO^L,div_free_fi)
    !w(ixO^S,nw+11)=div_free_fi
    !call get_Lorentz_force(ixI^L,ixO^L,w,Lorentz)
    !w(ixO^S,nw+12)=Lorentz(ixO^S,1)
    !w(ixO^S,nw+13)=Lorentz(ixO^S,2)
    !w(ixO^S,nw+14)=Lorentz(ixO^S,3)
    ! out dip
    !r (ixO^S) = x(ixO^S,1)
    !th(ixO^S) = x(ixO^S,2)
    !sinth(ixO^S) = sin(th(ixO^S))
    !eps_sin = 1.0d-6
    !where (abs(sinth(ixO^S)) < eps_sin) sinth(ixO^S) = eps_sin
    !
    !call gradient(w(ixI^S,mag(1)), ixI^L, ixO^L, 1, dBr_dr)
    !call gradient(w(ixI^S,mag(1)), ixI^L, ixO^L, 2, dBr_dth)
    !call gradient(w(ixI^S,mag(1)), ixI^L, ixO^L, 3, dBr_dph)
!
    !dip_dir(ixO^S) = w(ixO^S,mag(1))                        *dBr_dr(ixO^S)+&
    !                 w(ixO^S,mag(2))/r(ixO^S)               *dBr_dth(ixO^S)+&
    !                 w(ixO^S,mag(3))/(r(ixO^S)*sinth(ixO^S))*dBr_dr(ixO^S)
    !Bmag(ixO^S)  = sqrt(sum(w(ixO^S,mag(1:ndir))**2, dim=ndim+1))
    !epsBr(ixO^S) = 0.05d0*Bmag(ixO^S)
!
    !w(ixO^S,nw+13)=0.d0
    !where( abs(w(ixO^S,mag(1))) < epsBr(ixO^S) .and. dip_dir(ixO^S) >= 0.d0 )
    !  w(ixO^S,nw+13) = 1.d0
    !end where
!
    !w(ixO^S,nw+14) = dip_dir(ixO^S)
!
  end subroutine specialvar_output

  !==============================================================================
  ! Purpose: names for special variable output
  !==============================================================================
  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Br Bt Bp Jr Jt Jp divB pth Vr Lr Lt Lp'
  end subroutine specialvarnames_output

end module mod_usr
