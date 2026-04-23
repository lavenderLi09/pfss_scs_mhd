!> Analytical bipolar field in a Parker wind on a stretched spherical grid.
module mod_usr
  use mod_mhd
  implicit none
  double precision  :: k_B,miu0,mass_H,usr_grav,SRadius,rhob,Tiso
  double precision :: q_para,d_para,L_para
  double precision, allocatable :: ece4(:),ece3(:),ece2(:),ece1(:),ecz2(:),ecz1(:)
  double precision :: f_q, f_d, f_L
  double precision :: amr_r_core_max, amr_r_sheet_max, amr_theta_guard
  double precision :: amr_sheet_halfwidth, amr_br_ratio_max
  double precision :: amr_jb_eta_max_min, amr_jb_eta_mean_min
  double precision :: amr_joverb_base, amr_joverb_frac_refine, amr_joverb_frac_coarsen
  double precision :: amr_outer_relax_frac, amr_theta_relax_frac
  integer :: amr_force_max_level, amr_joverb_nhot_refine
  ! Parker solar wind parameters
  double precision  :: rc, Vs
  double precision  :: Vout, V_surface
  integer           :: ix^D
contains

  !==============================================================================
  ! Purpose: to include global parameters, set user methods, set coordinate 
  !          system and activate physics module.
  !==============================================================================
  subroutine usr_init()

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
    
    usr_grav=-2.74d4*unit_length/unit_velocity**2  ! solar gravity
    SRadius=6.955d10/unit_length                   ! Solar radius

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr
    usr_init_vector_potential=>initvecpot_usr
    usr_set_electric_field => boundary_electric_field
    usr_special_bc      => specialbound_usr
    usr_refine_grid     => special_refine_grid
    usr_aux_output     => specialvar_output
    usr_gravity         => gravity
    usr_add_aux_names  => specialvarnames_output
    usr_refine_threshold=> specialthreshold

    call set_coordinate_system('spherical_3D')
    call mhd_activate()

  end subroutine usr_init

  ! Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ f_q, f_d, f_L, &
         amr_r_core_max, amr_r_sheet_max, amr_theta_guard, &
         amr_sheet_halfwidth, amr_br_ratio_max, &
         amr_jb_eta_max_min, amr_jb_eta_mean_min, &
         amr_joverb_base, amr_joverb_frac_refine, amr_joverb_frac_coarsen, &
         amr_joverb_nhot_refine, &
         amr_force_max_level, &
         amr_outer_relax_frac, amr_theta_relax_frac

    f_q=1.d0
    f_d=1.d0
    f_L=1.d0
    amr_r_core_max=4.d0
    amr_r_sheet_max=8.d0
    amr_theta_guard=0.02d0
    amr_sheet_halfwidth=0.35d0
    amr_br_ratio_max=0.20d0
    amr_jb_eta_max_min=1.0d0
    amr_jb_eta_mean_min=0.3d0
    amr_joverb_base=0.02d0
    amr_joverb_frac_refine=0.05d0
    amr_joverb_frac_coarsen=0.01d0
    amr_joverb_nhot_refine=8
    amr_force_max_level=1
    amr_outer_relax_frac=0.25d0
    amr_theta_relax_frac=0.10d0
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  !==============================================================================
  ! Purpose: to initialize user public parameters and reset global parameters.
  !          Input data are also read here.
  !==============================================================================
  subroutine initglobaldata_usr()

    integer :: iv
    double precision :: qb,qc,qb23,qc23,q1b3,q1c3

    ! Set the analytical bipolar field.
    call usr_params_read(par_files)
    q_para= 1.0d20/(unit_magneticfield*unit_length**2)*f_q ! strength and sign ofmagnetic charges
    d_para= 1.0d9/unit_length*f_d ! depth of magnetic charges
    L_para= 2.0d9/unit_length*f_L ! half distance between magnetic charges
    if(mype==0) print*,'q_para',q_para*unit_magneticfield*unit_length**3
    if(mype==0) print*,'d_para',d_para*unit_length
    if(mype==0) print*,'L_para',L_para*unit_length
    
    ! Parker solar wind model.
    rc = 3.45d0
    Vs = 117.54d0 
    rhob = 5.0d9/1000.0d0/unit_numberdensity 
    Tiso = 2.0d0

    call cal_parker_solar_wind(1.0d0,rc,vs,V_surface)

    ! coefficent for gradient extrapolation in stretched mesh
    if(stretched_dim(1)) then
      allocate(ece4(refine_max_level),ece3(refine_max_level),ece2(refine_max_level),ece1(refine_max_level))
      allocate(ecz2(refine_max_level),ecz1(refine_max_level))
      do iv=1,refine_max_level
        qb=1.d0+qstretch(iv,1)
        qc=qb+qstretch(iv,1)**2
        qb23=qb**2-qb**3
        qc23=qc**2-qc**3
        q1b3=1.d0-qb**3
        q1c3=1.d0-qc**3
        ece1(iv)=1.d0/(qstretch(iv,1)*(q1b3*qc23-q1c3*qb23))
        ece4(iv)=qb23*ece1(iv)
        ece3(iv)=-(qstretch(iv,1)*qb23+qc23)*ece1(iv)
        ece2(iv)=(qb**3*qc23-qc**3*qb23+qstretch(iv,1)*qc23)*ece1(iv)
        ece1(iv)=(q1b3*qc23-q1c3*qb23-qstretch(iv,1)*(qb**3*qc23-qc**3*qb23))*ece1(iv)
        ecz2(iv)=1.d0/(1.d0-qb**2)
        ecz1(iv)=-qb**2*ecz2(iv)
      end do
    end if

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: A(ixI^S,1:ndim)
    double precision :: Bfr(ixI^S,1:ndim)
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
      print *,'Relax a solor wind model with politropic MHD'
      print *,'User Tiso(MK): ', Tiso
      print *,'V_surface(unit_velocity): ', V_surface
      first=.false.
      end if
    end if
    w(ixO^S,mom(1:3))=0.0d0
    if(stagger_grid) then ! I cant understand the CT part --------Yihua
      ! CT compute B by b_from_vector_potential, it use the subrutine `b_from_vector_potentialA` in mod_constrain_transport.t, and `b_from_vector_potential` call the user's subrutine `usr_init_vector_potential` is point to your `initvecpot_usr`. In this part, you can define the vector potential for the magnetic field you want. Then, the magnetic field compute by vector potential is in face due to the stagger_grid, to use FVM solver, you need to call `mhd_face_to_center` to transfer the quantity to the cell center. ----- Guoyin
      call b_from_vector_potential(block%ixGs^L,ixI^L,ixO^L,block%ws,x)
      call mhd_face_to_center(ixO^L,block)
    else
      call bipolar_field(ixI^L,ixO^L,x,A,Bfr)
      w(ixO^S,mag(:))=Bfr(ixO^S,:)
    end if

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
         call cal_parker_solar_wind(x(ix^D,1),rc,vs,Vout)
         w(ix^D,rho_)  = (rhob*V_surface)/(Vout*x(ix^D,1)**2)
         w(ix^D,mom(1))= Vout*w(ix^D,rho_)
    {end do\}
    w(ixO^S,p_) = w(ixO^S,rho_)*Tiso/(mhd_gamma-1.0d0) 

    !call mhd_to_conserved(ixI^L,ixO^L,w,x) 
  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    ! vector potential
    double precision :: Avec(ixI^S,1:ndim),Avec1(ixI^S,1:ndim), Avec2(ixI^S,1:ndim)
    ! Cartesian coordinates
    double precision :: xS(ixI^S,1:ndim)

    ! +++++++++ Special Regular :: 1=2, 2=3, 3=1 +++++++++
    xS(ixC^S,1) = xC(ixC^S,1) * cos(0.50d0*dpi-xC(ixC^S,2)) * cos(xC(ixC^S,3))
    xS(ixC^S,2) = xC(ixC^S,1) * cos(0.50d0*dpi-xC(ixC^S,2)) * sin(xC(ixC^S,3))
    xS(ixC^S,3) = xC(ixC^S,1) * sin(0.50d0*dpi-xC(ixC^S,2))
    Avec1=zero
    Avec2=zero
    call bipolar_field(ixI^L,ixC^L,xS,Avec1)
    Avec2=Avec1
    call Cart2SphereVector(ixI^L,ixC^L,xC,Avec2,Avec)

    if (idir==3) then
      A(ixC^S)=Avec(ixC^S,3)
    else if(idir==2) then 
      A(ixC^S)=Avec(ixC^S,2)
    else
      A(ixC^S)=Avec(ixC^S,1)
    end if

  end subroutine initvecpot_usr

  subroutine solve_sphere_curlvector(qvec,x,ixI^L,ixO^L,curlvec,&
                                     idirmin,idirmin0,ndir0)

    integer, intent(in)             :: ixI^L,ixO^L
    integer, intent(in)             :: ndir0, idirmin0
    integer, intent(inout)          :: idirmin
    double precision, intent(in)    :: qvec(ixI^S,1:ndir0)
    double precision, intent(in)    :: x(ixI^S,1:ndir0)
    double precision, intent(inout) :: curlvec(ixI^S,idirmin0:3)

    integer          :: ixA^L,ixC^L,jxC^L,idir,jdir,kdir,hxO^L,jxO^L,kxO^L,gxO^L
    double precision :: invdx(1:ndim)
    double precision :: tmp(ixI^S),tmp2(ixI^S),xC(ixI^S),surface(ixI^S)

    ! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
    ! Curl can have components (idirmin:3)
    ! Determine exact value of idirmin while doing the loop.
    ! Second order, stencil width is one
    ixA^L=ixO^L^LADD1;

    idirmin=4
    curlvec(ixI^S,idirmin0:3)=zero

    do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
      if(lvc(idir,jdir,kdir)/=0)then
        tmp(ixA^S)=qvec(ixA^S,kdir)
        hxO^L=ixO^L-kr(jdir,^D);
        jxO^L=ixO^L+kr(jdir,^D);
        select case(jdir)
          case(1)
            tmp(ixA^S)=tmp(ixA^S)*x(ixA^S,1)
            tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((x(jxO^S,1)-x(hxO^S,1))*x(ixO^S,1))
          case(2)
            if(idir==1) tmp(ixA^S)=tmp(ixA^S)*dsin(x(ixA^S,2))
            tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((x(jxO^S,2)-x(hxO^S,2))*x(ixO^S,1))
            if(idir==1) tmp2(ixO^S)=tmp2(ixO^S)/dsin(x(ixO^S,2))
          case(3)
            tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((x(jxO^S,3)-x(hxO^S,3))*&
                        x(ixO^S,1)*dsin(x(ixO^S,2)))
        end select
        if(lvc(idir,jdir,kdir)==1)then
          curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
        else
          curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
        endif
        if(idir<idirmin)idirmin=idir
      endif
    enddo; enddo; enddo;

  end subroutine solve_sphere_curlvector

  subroutine Cart2SphereVector(ixI^L,ixO^L,x,A_in,A_out)

    implicit none
    integer :: ixI^L,ixO^L
    double precision :: x(ixI^S,1:ndim)
    double precision :: A_in(ixI^S,1:ndim)
    double precision :: A_out(ixI^S,1:ndim)

    real*8 :: raddeg = dpi/ 180.
    real*8 :: lon(ixI^S),lat(ixI^S)
    real*8 :: bxCart(ixI^S),byCart(ixI^S),bzCart(ixI^S)
    real*8 :: br(ixI^S),bth(ixI^S),bph(ixI^S)
    real*8 :: a11(ixI^S),a12(ixI^S),a13(ixI^S)
    real*8 :: a21(ixI^S),a22(ixI^S),a23(ixI^S)
    real*8 :: a31(ixI^S),a32(ixI^S),a33(ixI^S)
    real*8 :: latc(ixI^S),lonc(ixI^S),pAng(ixI^S)

    bxCart=0.0d0
    byCart=0.0d0
    bzCart=0.0d0
    lon=0.0d0
    lat=0.0d0
    bph=0.0d0
    bth=0.0d0
     br=0.0d0
  A_out=0.0d0

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

  end subroutine Cart2SphereVector

  subroutine bipolar_field(ixI^L,ixO^L,x,A,Bbp)

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    ! vector potential
    double precision, intent(out)   :: A(ixI^S,1:ndim)
    ! magnetic field
    double precision, optional, intent(out)   :: Bbp(ixI^S,1:ndir)

    double precision :: Aphi(ixO^S),tmp(ixO^S)

    Aphi(ixO^S)= q_para*(L_para-x(ixO^S,3))/(sqrt(x(ixO^S,2)**2+(x(ixO^S,1)+d_para)**2)*&
                     sqrt(x(ixO^S,2)**2+(x(ixO^S,1)+d_para)**2+(x(ixO^S,3)-L_para)**2))+&
                 q_para*(L_para+x(ixO^S,3))/(sqrt(x(ixO^S,2)**2+(x(ixO^S,1)+d_para)**2)*&
                     sqrt(x(ixO^S,2)**2+(x(ixO^S,1)+d_para)**2+(x(ixO^S,3)+L_para)**2))

    A(ixO^S,1)=-Aphi(ixO^S)*x(ixO^S,2)/sqrt(x(ixO^S,2)**2+(x(ixO^S,1)+d_para)**2)
    A(ixO^S,3)= 0.d0
    A(ixO^S,2)= Aphi(ixO^S)*(x(ixO^S,1)+d_para)/sqrt(x(ixO^S,2)**2+(x(ixO^S,1)+d_para)**2)

    if(present(Bbp)) then
      tmp(ixO^S)=sqrt(x(ixO^S,2)**2+(x(ixO^S,1)+d_para)**2+(x(ixO^S,3)+L_para)**2)**3
      Bbp(ixO^S,3)=       (x(ixO^S,3)+L_para)/tmp(ixO^S)
      Bbp(ixO^S,2)=                x(ixO^S,2)/tmp(ixO^S)
      Bbp(ixO^S,1)= (x(ixO^S,1)+d_para)/tmp(ixO^S)
      tmp(ixO^S)=sqrt(x(ixO^S,2)**2+(x(ixO^S,1)+d_para)**2+(x(ixO^S,3)-L_para)**2)**3
      Bbp(ixO^S,3)=      -(x(ixO^S,3)-L_para)/tmp(ixO^S) + Bbp(ixO^S,3)
      Bbp(ixO^S,2)=               -x(ixO^S,2)/tmp(ixO^S) + Bbp(ixO^S,2)
      Bbp(ixO^S,1)=-(x(ixO^S,1)+d_para)/tmp(ixO^S) + Bbp(ixO^S,1)
      Bbp(ixO^S,:)=q_para*Bbp(ixO^S,:)
    end if

  end subroutine bipolar_field

  subroutine set_analytic_parker_bipole(ixI^L,ixO^L,x,w)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: A(ixI^S,1:ndim),Bcart(ixI^S,1:ndim),Bsph(ixI^S,1:ndim),xS(ixI^S,1:ndim)

    w(ixO^S,mom(:))=0.d0
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      call cal_parker_solar_wind(x(ix^D,1), rc, vs, vout)
      w(ix^D,rho_)  = (rhob*V_surface)/(Vout*x(ix^D,1)**2)
      w(ix^D,mom(1))= vout*w(ix^D,rho_)
      w(ix^D,p_)    = w(ix^D,rho_)*Tiso/(mhd_gamma-1.0d0)
    {end do\}

    xS(ixO^S,1) = x(ixO^S,1) * cos(0.50d0*dpi-x(ixO^S,2)) * cos(x(ixO^S,3))
    xS(ixO^S,2) = x(ixO^S,1) * cos(0.50d0*dpi-x(ixO^S,2)) * sin(x(ixO^S,3))
    xS(ixO^S,3) = x(ixO^S,1) * sin(0.50d0*dpi-x(ixO^S,2))
    call bipolar_field(ixI^L,ixO^L,xS,A,Bcart)
    call Cart2SphereVector(ixI^L,ixO^L,x,Bcart,Bsph)
    w(ixO^S,mag(:))=Bsph(ixO^S,:)
  end subroutine set_analytic_parker_bipole

  subroutine boundary_electric_field(ixI^L,ixO^L,qt,qdt,fE,s)
    ! specify tangential electric field at physical boundaries 
    ! to fix or drive normal magnetic field
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt,qdt
    type(state)                        :: s
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    double precision :: xC(ixI^S,1:ndim),xlen^D
    integer :: idir,ixC^L,ixA^L,ix^D,ixg^D

    associate(w=>s%w,x=>s%x,ws=>s%ws)

    if(s%is_physical_boundary(1)) then
      ixCmin^D=ixOmin^D-1;
      ixCmax^D=ixOmax^D;
      fE(nghostcells^%1ixC^S,2:3)=0.0d0
    end if

    ! perfect conductor theta boundaries
    ! is this needed if we want a pole boundary at theta? -------Yihua
    if(block%is_physical_boundary(3)) then
      ixCmin^D=ixOmin^D-kr(2,^D);
      ixCmax^D=ixOmax^D;
      ixAmin^D=ixCmin^D-kr(3,^D);
      ixAmax^D=ixCmax^D;
      fE(ixAmin2^%2ixA^S,1)=0.d0
      ixAmin^D=ixCmin^D-kr(1,^D);
      ixAmax^D=ixCmax^D;
      fE(ixAmin2^%2ixA^S,3)=0.d0
    end if
    if(block%is_physical_boundary(4)) then
      ixCmin^D=ixOmin^D-kr(2,^D);
      ixCmax^D=ixOmax^D;
      ixAmin^D=ixCmin^D-kr(3,^D);
      ixAmax^D=ixCmax^D;
      fE(ixAmax2^%2ixA^S,1)=0.d0
      ixAmin^D=ixCmin^D-kr(1,^D);
      ixAmax^D=ixCmax^D;
      fE(ixAmax2^%2ixA^S,3)=0.d0
    end if

    end associate

  end subroutine boundary_electric_field

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

  end subroutine cal_parker_solar_wind

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    
    double precision :: Qp(ixI^S)
    double precision :: A(ixI^S,1:ndim),Bcart(ixI^S,1:ndim),Bsph(ixI^S,1:ndim),xS(ixI^S,1:ndim)
    integer :: ix^D,ixOs^L,jxO^L,idir
    double precision :: q,q1,b0,b1,b2
    
    q  = dble(qstretch_baselevel(1))
    q1 = 1.0d0/q
    b0 = -(2.0d0+q1)/(1.0d0+q1)
    b1 = 1.0d0+1.0d0/q1
    b2 = -1.0d0/(1.0d0+q1**2)

    select case(iB)

 case(1)
       w(ixO^S,mom(:))=0.d0
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
         call cal_parker_solar_wind(x(ix^D,1), rc, vs, vout)
         w(ix^D,rho_)  = (rhob*V_surface)/(Vout*x(ix^D,1)**2)
         w(ix^D,mom(1))= vout*w(ix^D,rho_)
         w(ix^D,p_)    = w(ix^D,rho_)*Tiso/(mhd_gamma-1.0d0)
       {end do\}

       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           if(stretched_dim(1)) then
             do ix1=ixOsmax1,ixOsmin1,-1
                block%ws(ix1^%1ixOs^S,idir) = &
                   ece4(block%level)*block%ws(ix1+4^%1ixOs^S,idir)+&
                   ece3(block%level)*block%ws(ix1+3^%1ixOs^S,idir)+&
                   ece2(block%level)*block%ws(ix1+2^%1ixOs^S,idir)+&
                   ece1(block%level)*block%ws(ix1+1^%1ixOs^S,idir)
             end do
           else 
             do ix1=ixOsmax1,ixOsmin1,-1
              block%ws(ix1^%1ixOs^S,idir) = &
                  0.12d0*block%ws(ix1+5^%1ixOs^S,idir)&
                 -0.76d0*block%ws(ix1+4^%1ixOs^S,idir)&
                 +2.08d0*block%ws(ix1+3^%1ixOs^S,idir)&
                 -3.36d0*block%ws(ix1+2^%1ixOs^S,idir)&
                 +2.92d0*block%ws(ix1+1^%1ixOs^S,idir)

              !block%ws(ix1^%1ixOs^S,idir) = third*&
              !       (-block%ws(ix1+2^%1ixOs^S,idir)&
              !   +4.d0*block%ws(ix1+1^%1ixOs^S,idir))
             end do
            end if 
         end do
         ixOs^L=ixO^L-kr(1,^D);
         jxO^L=ixO^L+nghostcells*kr(1,^D);
         block%ws(ixOs^S,1)=zero
         do ix1=ixOsmax1,ixOsmin1,-1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix1^%1ixOs^S,1)=Qp(ix1+1^%1ixO^S)*block%dvolume(ix1+1^%1ixO^S)&
             /block%surfaceC(ix1^%1ixOs^S,1)
         end do
         call mhd_face_to_center(ixO^L,block)

         ! Keep the cell-centered ghost magnetic field consistent with the
         ! analytical buried-charge bipole used in the interior initialization.
         xS(ixO^S,1) = x(ixO^S,1) * cos(0.50d0*dpi-x(ixO^S,2)) * cos(x(ixO^S,3))
         xS(ixO^S,2) = x(ixO^S,1) * cos(0.50d0*dpi-x(ixO^S,2)) * sin(x(ixO^S,3))
         xS(ixO^S,3) = x(ixO^S,1) * sin(0.50d0*dpi-x(ixO^S,2))
         call bipolar_field(ixI^L,ixO^L,xS,A,Bcart)
         call Cart2SphereVector(ixI^L,ixO^L,x,Bcart,Bsph)
         w(ixO^S,mag(:))=Bsph(ixO^S,:)
       else
         do ix1=ixOmax1,ixOmin1,-1
           w(ix1^%1ixO^S,mag(:))=third* &
                      (-w(ix1+2^%1ixO^S,mag(:)) &
                 +4.0d0*w(ix1+1^%1ixO^S,mag(:)))
         enddo
        end if

 case(2)
       !do ix1=ixOmin1,ixOmax1
       !  w(ix1^%1ixO^S,mom(1))=w(ixOmin1-1^%1ixO^S,mom(1))
       !  w(ix1^%1ixO^S,mom(2))=w(ixOmin1-1^%1ixO^S,mom(2))
       !  w(ix1^%1ixO^S,mom(3))=w(ixOmin1-1^%1ixO^S,mom(3))
       !enddo
       !w(ixO^S,mom(1))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))
       !w(ixO^S,mom(2))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))
       !w(ixO^S,mom(3))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))
       
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
         end do

       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
             ixOsmax^D=ixOmax^D;
             ixOsmin^D=ixOmin^D-kr(^D,idir);
             do ix1=ixOsmin1,ixOsmax1
                block%ws(ix1^%1ixOs^S,idir) = 0.0d0
                !block%ws(ix1^%1ixOs^S,idir) = 1.d0/11.d0*&
                !   ( -2.d0*block%ws(ix1-4^%1ixOs^S,idir)&
                !    +11.d0*block%ws(ix1-3^%1ixOs^S,idir)&
                !    -27.d0*block%ws(ix1-2^%1ixOs^S,idir)&
                !    +29.d0*block%ws(ix1-1^%1ixOs^S,idir))
                !block%ws(ix1^%1ixOs^S,idir) = third*&
                !       (-block%ws(ix1-2^%1ixOs^S,idir)&
                !   +4.d0*block%ws(ix1-1^%1ixOs^S,idir))
             end do
         end do
         ixOs^L=ixO^L;
         jxO^L=ixO^L-nghostcells*kr(1,^D);
         block%ws(ixOs^S,1)=zero
         do ix1=ixOsmin1,ixOsmax1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix1^%1ixOs^S,1)=-Qp(ix1^%1ixO^S)*block%dvolume(ix1^%1ixO^S)&
             /block%surfaceC(ix1^%1ixOs^S,1)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         do ix1=ixOmin1,ixOmax1
           w(ix1^%1ixO^S,mag(1))=third* &
                      (-w(ix1-2^%1ixO^S,mag(1)) &
                 +4.0d0*w(ix1-1^%1ixO^S,mag(1)))
         enddo
         w(ixO^S,mag(2:3))=0.d0
       end if

     case(3)
!      do ix2=ixOmax2,ixOmin2,-1
!        w(ix2^%2ixO^S,rho_)   = w(ix2+1^%2ixO^S,rho_)
!        w(ix2^%2ixO^S,p_)     = w(ix2+1^%2ixO^S,p_)
!        w(ix2^%2ixO^S,mom(:)) = w(ix2+1^%2ixO^S,mom(:))
!      end do
!
!      if(stagger_grid) then
!        do idir=1,nws
!          if(idir==2) cycle
!          ixOsmax^D=ixOmax^D;
!          ixOsmin^D=ixOmin^D-kr(^D,idir);
!          do ix2=ixOsmax2,ixOsmin2,-1
!            block%ws(ix2^%2ixOs^S,idir)=third*&
!                   (-block%ws(ix2+2^%2ixOs^S,idir)&
!               +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
!          end do
!        end do
!        ixOs^L=ixO^L-kr(2,^D);
!        jxO^L=ixO^L+nghostcells*kr(2,^D);
!        block%ws(ixOs^S,2)=zero
!        do ix2=ixOsmax2,ixOsmin2,-1
!          call get_divb(w,ixI^L,ixO^L,Qp)
!          block%ws(ix2^%2ixOs^S,2)=Qp(ix2+1^%2ixO^S)*block%dvolume(ix2+1^%2ixO^S)&
!            /block%surfaceC(ix2^%2ixOs^S,2)
!        end do
!        call mhd_face_to_center(ixO^L,block)
!
!        xS(ixO^S,1) = x(ixO^S,1) * cos(0.50d0*dpi-x(ixO^S,2)) * cos(x(ixO^S,3))
!        xS(ixO^S,2) = x(ixO^S,1) * cos(0.50d0*dpi-x(ixO^S,2)) * sin(x(ixO^S,3))
!        xS(ixO^S,3) = x(ixO^S,1) * sin(0.50d0*dpi-x(ixO^S,2))
!        call bipolar_field(ixI^L,ixO^L,xS,A,Bcart)
!        call Cart2SphereVector(ixI^L,ixO^L,x,Bcart,Bsph)
!        w(ixO^S,mag(:))=Bsph(ixO^S,:)
!      else
!        do ix2=ixOmax2,ixOmin2,-1
!          w(ix2^%2ixO^S,mag(:))=third*&
!                 (-w(ix2+2^%2ixO^S,mag(:))&
!             +4.d0*w(ix2+1^%2ixO^S,mag(:)))
!        end do
!      end if

       ! Analytic continuation across the truncated theta boundary:
       ! extend the same Parker wind and analytical bipolar field into the
       ! ghost cells instead of mixing zero-gradient hydro with extrapolated B.
       call set_analytic_parker_bipole(ixI^L,ixO^L,x,w)

       if(stagger_grid) then
         ! The centered analytical overwrite is fine for ghost cells, but the
         ! CT theta-cut faces still need the divergence-consistent reconstruction.
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix2=ixOsmax2,ixOsmin2,-1
             block%ws(ix2^%2ixOs^S,idir)=third*&
                    (-block%ws(ix2+2^%2ixOs^S,idir)&
                +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L-kr(2,^D);
         jxO^L=ixO^L+nghostcells*kr(2,^D);
         block%ws(ixOs^S,2)=zero
         do ix2=ixOsmax2,ixOsmin2,-1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix2^%2ixOs^S,2)=Qp(ix2+1^%2ixO^S)*block%dvolume(ix2+1^%2ixO^S)&
             /block%surfaceC(ix2^%2ixOs^S,2)
         end do
         call mhd_face_to_center(ixO^L,block)
         call set_analytic_parker_bipole(ixI^L,ixO^L,x,w)
       end if

     case(4)
!      do ix2=ixOmin2,ixOmax2
!        w(ix2^%2ixO^S,rho_)   = w(ix2-1^%2ixO^S,rho_)
!        w(ix2^%2ixO^S,p_)     = w(ix2-1^%2ixO^S,p_)
!        w(ix2^%2ixO^S,mom(:)) = w(ix2-1^%2ixO^S,mom(:))
!      end do
!
!      if(stagger_grid) then
!        do idir=1,nws
!          if(idir==2) cycle
!          ixOsmax^D=ixOmax^D;
!          ixOsmin^D=ixOmin^D-kr(^D,idir);
!          do ix2=ixOsmin2,ixOsmax2
!            block%ws(ix2^%2ixOs^S,idir)=third*&
!                   (-block%ws(ix2-2^%2ixOs^S,idir)&
!               +4.d0*block%ws(ix2-1^%2ixOs^S,idir))
!          end do
!        end do
!        ixOs^L=ixO^L;
!        jxO^L=ixO^L-nghostcells*kr(2,^D);
!        block%ws(ixOs^S,2)=zero
!        do ix2=ixOsmin2,ixOsmax2
!          call get_divb(w,ixI^L,ixO^L,Qp)
!          block%ws(ix2^%2ixOs^S,2)=-Qp(ix2^%2ixO^S)*block%dvolume(ix2^%2ixO^S)&
!            /block%surfaceC(ix2^%2ixOs^S,2)
!        end do
!        call mhd_face_to_center(ixO^L,block)
!
!        xS(ixO^S,1) = x(ixO^S,1) * cos(0.50d0*dpi-x(ixO^S,2)) * cos(x(ixO^S,3))
!        xS(ixO^S,2) = x(ixO^S,1) * cos(0.50d0*dpi-x(ixO^S,2)) * sin(x(ixO^S,3))
!        xS(ixO^S,3) = x(ixO^S,1) * sin(0.50d0*dpi-x(ixO^S,2))
!        call bipolar_field(ixI^L,ixO^L,xS,A,Bcart)
!        call Cart2SphereVector(ixI^L,ixO^L,x,Bcart,Bsph)
!        w(ixO^S,mag(:))=Bsph(ixO^S,:)
!      else
!        do ix2=ixOmin2,ixOmax2
!          w(ix2^%2ixO^S,mag(:))=third*&
!                 (-w(ix2-2^%2ixO^S,mag(:))&
!             +4.d0*w(ix2-1^%2ixO^S,mag(:)))
!        end do
!      end if

       ! Analytic continuation across the truncated theta boundary:
       ! extend the same Parker wind and analytical bipolar field into the
       ! ghost cells instead of mixing zero-gradient hydro with extrapolated B.
       call set_analytic_parker_bipole(ixI^L,ixO^L,x,w)

       if(stagger_grid) then
         ! The centered analytical overwrite is fine for ghost cells, but the
         ! CT theta-cut faces still need the divergence-consistent reconstruction.
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix2=ixOsmin2,ixOsmax2
             block%ws(ix2^%2ixOs^S,idir)=third*&
                    (-block%ws(ix2-2^%2ixOs^S,idir)&
                +4.d0*block%ws(ix2-1^%2ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L;
         jxO^L=ixO^L-nghostcells*kr(2,^D);
         block%ws(ixOs^S,2)=zero
         do ix2=ixOsmin2,ixOsmax2
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix2^%2ixOs^S,2)=-Qp(ix2^%2ixO^S)*block%dvolume(ix2^%2ixO^S)&
             /block%surfaceC(ix2^%2ixOs^S,2)
         end do
         call mhd_face_to_center(ixO^L,block)
         call set_analytic_parker_bipole(ixI^L,ixO^L,x,w)
       end if
     case(5)
   
     ! phi boundaries are defined by the periodic condition in AMRVAC

     case(6)
  
     ! phi boundaries are defined by the periodic condition in AMRVAC

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
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision :: rmin, rmax
    double precision :: dr_block, dtheta_block, dphi_block, dxmin_phys
    double precision :: joverb_cr, frac_hot
    double precision :: bmag(ixI^S), current_mag(ixI^S), delta_min(ixI^S), joverb(ixI^S)
    double precision :: current(ixI^S,3)
    logical :: dbg_print
    integer :: idirmin, nhot, ntheta, target_level

    refine=-1
    coarsen=-1

    rmin=minval(x(ixO^S,1))
    rmax=maxval(x(ixO^S,1))
    if(rmin >= amr_r_sheet_max) return

    ntheta=count(abs(x(ixO^S,2)-0.5d0*(xprobmin2+xprobmax2)) <= amr_sheet_halfwidth)
    if(ntheta <= 0) return

    dr_block=(maxval(x(ixO^S,1))-minval(x(ixO^S,1)))/dble(max(1,ixOmax1-ixOmin1))
    dtheta_block=(maxval(x(ixO^S,2))-minval(x(ixO^S,2)))/dble(max(1,ixOmax2-ixOmin2))
    dphi_block=(maxval(x(ixO^S,3))-minval(x(ixO^S,3)))/dble(max(1,ixOmax3-ixOmin3))

    dxlevel(1)=rnode(rpdx1_,igrid)
    dxlevel(2)=rnode(rpdx2_,igrid)
    dxlevel(3)=rnode(rpdx3_,igrid)
    block=>ps(igrid)

    ! get_current only guarantees valid values on ixO, so keep derived AMR
    ! diagnostics on the same active stencil instead of touching ghost cells.
    bmag(ixO^S)=sqrt(w(ixO^S,mag(1))**2 + w(ixO^S,mag(2))**2 + w(ixO^S,mag(3))**2)
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    current_mag(ixO^S)=sqrt(current(ixO^S,1)**2 + current(ixO^S,2)**2 + current(ixO^S,3)**2)
    delta_min(ixO^S)=min(dr_block, x(ixO^S,1)*dtheta_block, &
         x(ixO^S,1)*max(dsin(x(ixO^S,2)),smalldouble)*dphi_block)
    joverb(ixO^S)=current_mag(ixO^S)/(bmag(ixO^S)+smalldouble)

    dxmin_phys=max(minval(delta_min(ixO^S)), smalldouble)
    joverb_cr=amr_joverb_base/dxmin_phys
    nhot=count((abs(x(ixO^S,2)-0.5d0*(xprobmin2+xprobmax2)) <= amr_sheet_halfwidth) .and. &
               (joverb(ixO^S) >= joverb_cr))
    frac_hot=dble(nhot)/dble(ntheta)

    if(rmax <= amr_r_core_max) then
      target_level=min(refine_max_level,2)
    else
      target_level=min(refine_max_level,1)
    end if

    dbg_print = (nhot > 0) .or. (level /= target_level) .or. &
                (frac_hot >= 0.2d0*amr_joverb_frac_refine)

    if(level > target_level) then
      coarsen=1
      if(dbg_print) then
        write(*,*) 'AMRDBG1',mype,igrid,level,rmin,rmax,joverb_cr,frac_hot
        write(*,*) 'AMRDBG2',ntheta,nhot,target_level,refine,coarsen
      end if
      return
    end if

    if(level < target_level) then
      if(frac_hot >= amr_joverb_frac_refine .or. nhot >= amr_joverb_nhot_refine) then
        refine=1
        coarsen=-1
      end if
    end if

    if(refine == 1) then
      if(dbg_print) then
        write(*,*) 'AMRDBG1',mype,igrid,level,rmin,rmax,joverb_cr,frac_hot
        write(*,*) 'AMRDBG2',ntheta,nhot,target_level,refine,coarsen
      end if
      return
    end if

    if(frac_hot <= amr_joverb_frac_coarsen) then
      if(level > 1) then
        coarsen=1
      end if
      if(dbg_print) then
        write(*,*) 'AMRDBG1',mype,igrid,level,rmin,rmax,joverb_cr,frac_hot
        write(*,*) 'AMRDBG2',ntheta,nhot,target_level,refine,coarsen
      end if
      return
    end if

    if(dbg_print) then
      write(*,*) 'AMRDBG1',mype,igrid,level,rmin,rmax,joverb_cr,frac_hot
      write(*,*) 'AMRDBG2',ntheta,nhot,target_level,refine,coarsen
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

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: ws1(ixGs^T,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: w_om(ixI^S,nw)
    double precision                   :: qvec(ixI^S,1:ndim)
    double precision                   :: bvec(ixI^S,1:ndim)
    double precision                   :: patchwi(ixI^S)
    double precision                   :: xO^L
    double precision                   :: tmp(ixI^S),tmp1(ixI^S)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: csound2(ixI^S)

    double precision :: current(ixI^S,3)
    double precision :: xlen^D
    integer :: ix^D,ixbc^D,idirmin,idir,jdir,kdir,hxO^L,idim

    w(ixO^S,nw+1)=w(ixO^S,mag(1))
    w(ixO^S,nw+2)=w(ixO^S,mag(2))
    w(ixO^S,nw+3)=w(ixO^S,mag(3))
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+4)=current(ixO^S,1)
    w(ixO^S,nw+5)=current(ixO^S,2)
    w(ixO^S,nw+6)=current(ixO^S,3)
    {xOmin^D = xprobmin^D + 0.05d0*(xprobmax^D-xprobmin^D)\}
    {xOmax^D = xprobmax^D - 0.05d0*(xprobmax^D-xprobmin^D)\}
    xOmin1 = xprobmin1

    qvec(ixI^S,1:ndir)=zero
    bvec(ixI^S,1:ndir)=zero
    tmp(ixI^S)=zero
    tmp1(ixI^S)=zero
    w_om=zero
    patchwi=zero

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if({ x(ix^DD,^D) > xOmin^D .and. x(ix^DD,^D) < xOmax^D | .and. }) then
          patchwi(ix^D)=1.0d0
        else
          patchwi(ix^D)=0.0d0
        endif
    {end do\}

    if(B0field) then
      bvec(ixI^S,:)=w(ixI^S,mag(:))+block%b0(ixI^S,mag(:),0)
    else
      bvec(ixI^S,:)=w(ixI^S,mag(:))
    endif
    do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
       if(lvc(idir,jdir,kdir)/=0)then
          tmp(ixO^S)=current(ixO^S,jdir)*bvec(ixO^S,kdir)
          if(lvc(idir,jdir,kdir)==1)then
             qvec(ixO^S,idir)=qvec(ixO^S,idir)+tmp(ixO^S)
          else
             qvec(ixO^S,idir)=qvec(ixO^S,idir)-tmp(ixO^S)
          endif
       endif
    enddo; enddo; enddo

    w(ixO^S,nw+7)=patchwi(ixO^S)*sqrt(sum(qvec(ixO^S,:)**2,ndim+1))/&
                   sqrt(sum(bvec(ixO^S,:)**2,ndim+1))
    w_om(ixO^S,1:nw)=w(ixO^S,1:nw)
    call get_normalized_divb(w_om,ixI^L,ixO^L,tmp1)
    w(ixO^S,nw+8)=tmp1(ixO^S) 
    w(ixO^S,nw+9)=w(ixO^S,mom(1))*unit_velocity/(w(ixO^S,rho_)*1.0d5)
    w(ixO^S,nw+10)=w(ixO^S,mom(2))*unit_velocity/(w(ixO^S,rho_)*1.0d5)
    w(ixO^S,nw+11)=w(ixO^S,mom(3))*unit_velocity/(w(ixO^S,rho_)*1.0d5)

  end subroutine specialvar_output

  !==============================================================================
  ! Purpose: names for special variable output
  !==============================================================================
  subroutine specialvarnames_output(varnames)
    use mod_global_parameters

    character(len=*) varnames

    varnames='br bt bp J1 J2 J3 CW fi Vr Vt Vp'

  end subroutine specialvarnames_output

  subroutine specialthreshold(wlocal,xlocal,tolerance,qt,level)
    !PURPOSE: use different tolerance in special regions for AMR to
    !reduce/increase resolution there where nothing/something interesting happens.
    use mod_global_parameters

    double precision, intent(in) :: wlocal(1:nw),xlocal(1:ndim),qt
    double precision, intent(inout) :: tolerance
    integer, intent(in) :: level
    
    double precision :: outer_zone, tol_add

    tol_add=0.d0
    outer_zone=amr_outer_relax_frac*(xprobmax1-xprobmin1)

    if(xprobmax1-xlocal(1) < outer_zone) then
      tol_add=tol_add + (1.d0-(xprobmax1-xlocal(1))/outer_zone)*0.5d0
    endif

    tolerance=tolerance+tol_add

  end subroutine specialthreshold

end module mod_usr
