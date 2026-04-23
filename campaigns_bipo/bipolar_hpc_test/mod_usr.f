!> Analytical bipolar field in a Parker wind on a stretched spherical grid.
module mod_usr
  use mod_mhd
  implicit none
  double precision  :: k_B,miu0,mass_H,usr_grav,SRadius,rhob,Tiso
  double precision :: q_para,d_para,L_para
  double precision, allocatable :: ece4(:),ece3(:),ece2(:),ece1(:),ecz2(:),&
     ecz1(:)
  double precision :: f_q, f_d, f_L
  double precision :: amr_r_core_max, amr_r_sheet_max, amr_theta_guard
  double precision :: amr_sheet_halfwidth, amr_br_ratio_max
  double precision :: amr_jb_eta_max_min, amr_jb_eta_mean_min
  double precision :: amr_outer_relax_frac, amr_theta_relax_frac
  integer :: amr_force_max_level
  ! Parker solar wind parameters
  double precision  :: rc, Vs
  double precision  :: Vout, V_surface
  integer           :: ix1,ix2,ix3
contains

  !==============================================================================
  ! Purpose: to include global parameters, set user methods, set coordinate 
  !          system and activate physics module.
  !==============================================================================
  subroutine usr_init()

    ! normalization unit in CGS Unit
    k_B = 1.3806d-16          ! erg*K-1,erg*K-1,erg*K-1
    miu0 = 4.d0*dpi !Gauss2,Gauss2,Gauss2 cm2,cm2,cm2 dyne-1,dyne-1,dyne-1
    mass_H = 1.67262d-24      ! g
    unit_length        = 6.955d10 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm-3,cm-3,cm-3
    unit_density       = 1.4d0*mass_H*unit_numberdensity !2.341668000000000E-015 g*cm-3,g*cm-3,g*cm-3
    unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature !0.317538000000000 erg*cm-3,erg*cm-3,erg*cm-3
    unit_magneticfield = dsqrt(miu0*unit_pressure) !1.99757357615242 Gauss
    unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density) !1.16448846777562E007 cm/s = 116.45 km/s
    unit_time          = unit_length/unit_velocity !5972.5794 s = 99.543 min
    
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

    namelist /usr_list/ f_q, f_d, f_L, amr_r_core_max, amr_r_sheet_max,&
        amr_theta_guard, amr_sheet_halfwidth, amr_br_ratio_max,&
        amr_jb_eta_max_min, amr_jb_eta_mean_min, amr_force_max_level,&
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
    q_para= 1.0d20/(unit_magneticfield*unit_length**2)*f_q !strength and sign ofmagnetic charges
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
      allocate(ece4(refine_max_level),ece3(refine_max_level),&
         ece2(refine_max_level),ece1(refine_max_level))
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
        ece1(iv)=(q1b3*qc23-q1c3*qb23-qstretch(iv,&
           1)*(qb**3*qc23-qc**3*qb23))*ece1(iv)
        ecz2(iv)=1.d0/(1.d0-qb**2)
        ecz1(iv)=-qb**2*ecz2(iv)
      end do
    end if

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: A(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: Bfr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
      print *,'Relax a solor wind model with politropic MHD'
      print *,'User Tiso(MK): ', Tiso
      print *,'V_surface(unit_velocity): ', V_surface
      first=.false.
      end if
    end if
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1:3))=0.0d0
    if(stagger_grid) then ! I cant understand the CT part --------Yihua
      ! CT compute B by b_from_vector_potential, it use the subrutine `b_from_vector_potentialA` in mod_constrain_transport.t, and `b_from_vector_potential` call the user's subrutine `usr_init_vector_potential` is point to your `initvecpot_usr`. In this part, you can define the vector potential for the magnetic field you want. Then, the magnetic field compute by vector potential is in face due to the stagger_grid, to use FVM solver, you need to call `mhd_face_to_center` to transfer the quantity to the cell center. ----- Guoyin
      call b_from_vector_potential(block%ixGsmin1,block%ixGsmin2,&
         block%ixGsmin3,block%ixGsmax1,block%ixGsmax2,block%ixGsmax3,ixImin1,&
         ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
         ixOmax1,ixOmax2,ixOmax3,block%ws,x)
      call mhd_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         block)
    else
      call bipolar_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,A,Bfr)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))=Bfr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)
    end if

    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
         call cal_parker_solar_wind(x(ix1,ix2,ix3,1),rc,vs,Vout)
         w(ix1,ix2,ix3,rho_)  = (rhob*V_surface)/(Vout*x(ix1,ix2,ix3,1)**2)
         w(ix1,ix2,ix3,mom(1))= Vout*w(ix1,ix2,ix3,rho_)
    end do
    end do
    end do
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)*Tiso/(mhd_gamma-1.0d0) 

    !call mhd_to_conserved(ixI^L,ixO^L,w,x) 
  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idir
    double precision, intent(in)       :: xC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)      :: A(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    ! vector potential
    double precision :: Avec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim),Avec1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim),&
        Avec2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    ! Cartesian coordinates
    double precision :: xS(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)

    ! +++++++++ Special Regular :: 1=2, 2=3, 3=1 +++++++++
    xS(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1) = xC(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,ixCmin3:ixCmax3,1) * cos(0.50d0*dpi-xC(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)) * cos(xC(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,ixCmin3:ixCmax3,3))
    xS(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2) = xC(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,ixCmin3:ixCmax3,1) * cos(0.50d0*dpi-xC(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)) * sin(xC(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,ixCmin3:ixCmax3,3))
    xS(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,3) = xC(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,ixCmin3:ixCmax3,1) * sin(0.50d0*dpi-xC(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,ixCmin3:ixCmax3,2))
    Avec1=zero
    Avec2=zero
    call bipolar_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
       ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,xS,Avec1)
    Avec2=Avec1
    call Cart2SphereVector(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,xC,Avec2,Avec)

    if (idir==3) then
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=Avec(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,3)
    else if(idir==2) then 
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=Avec(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)
    else
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=Avec(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)
    end if

  end subroutine initvecpot_usr

  subroutine solve_sphere_curlvector(qvec,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,curlvec,&
     idirmin,idirmin0,ndir0)

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    integer, intent(in)             :: ndir0, idirmin0
    integer, intent(inout)          :: idirmin
    double precision, intent(in)    :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir0)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir0)
    double precision, intent(inout) :: curlvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,idirmin0:3)

    integer          :: ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,&
       ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,jxCmin1,jxCmin2,jxCmin3,&
       jxCmax1,jxCmax2,jxCmax3,idir,jdir,kdir,hxOmin1,hxOmin2,hxOmin3,hxOmax1,&
       hxOmax2,hxOmax3,jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3,kxOmin1,&
       kxOmin2,kxOmin3,kxOmax1,kxOmax2,kxOmax3,gxOmin1,gxOmin2,gxOmin3,gxOmax1,&
       gxOmax2,gxOmax3
    double precision :: invdx(1:ndim)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       xC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       surface(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    ! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
    ! Curl can have components (idirmin:3)
    ! Determine exact value of idirmin while doing the loop.
    ! Second order, stencil width is one
    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;

    idirmin=4
    curlvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idirmin0:3)=zero

    do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
      if(lvc(idir,jdir,kdir)/=0)then
        tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3)=qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,kdir)
        hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
        hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
        hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
        jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
        jxOmin3=ixOmin3+kr(jdir,3);jxOmax1=ixOmax1+kr(jdir,1)
        jxOmax2=ixOmax2+kr(jdir,2);jxOmax3=ixOmax3+kr(jdir,3);
        select case(jdir)
          case(1)
            tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3)=tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3)*x(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3,1)
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
               jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3))/((x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
               jxOmin3:jxOmax3,1)-x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3,1))*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,1))
          case(2)
            if(idir==1) tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3)=tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3)*dsin(x(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3,2))
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
               jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3))/((x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
               jxOmin3:jxOmax3,2)-x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3,2))*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,1))
            if(idir==1) tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,2))
          case(3)
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
               jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3))/((x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
               jxOmin3:jxOmax3,3)-x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3,3))*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,2)))
        end select
        if(lvc(idir,jdir,kdir)==1)then
          curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        else
          curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        endif
        if(idir<idirmin)idirmin=idir
      endif
    enddo; enddo; enddo;

  end subroutine solve_sphere_curlvector

  subroutine Cart2SphereVector(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,A_in,A_out)

    implicit none
    integer :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision :: x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: A_in(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: A_out(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)

    real*8 :: raddeg = dpi/ 180.
    real*8 :: lon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       lat(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real*8 :: bxCart(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       byCart(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       bzCart(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real*8 :: br(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       bth(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       bph(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real*8 :: a11(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       a12(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       a13(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real*8 :: a21(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       a22(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       a23(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real*8 :: a31(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       a32(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       a33(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real*8 :: latc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       lonc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       pAng(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    bxCart=0.0d0
    byCart=0.0d0
    bzCart=0.0d0
    lon=0.0d0
    lat=0.0d0
    bph=0.0d0
    bth=0.0d0
     br=0.0d0
  A_out=0.0d0

    bxCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = A_in(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)
    byCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = A_in(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)
    bzCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = A_in(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)
    lon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) =  x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
    lat(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) =  0.5d0*dpi - &
       x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
    latc = 0.0d0
    lonc = 0.0d0
    pAng = 0.0d0

    a11 = -sin(latc) * sin(pAng) * sin(lon - lonc) + cos(pAng) * cos(lon - &
       lonc)
    a12 =  sin(latc) * cos(pAng) * sin(lon - lonc) + sin(pAng) * cos(lon - &
       lonc)
    a13 = -cos(latc) * sin(lon - lonc)
    a21 = -sin(lat) * (sin(latc) * sin(pAng) * cos(lon - lonc) + cos(pAng) * &
       sin(lon - lonc)) - cos(lat) * cos(latc) * sin(pAng)
    a22 =  sin(lat) * (sin(latc) * cos(pAng) * cos(lon - lonc) - sin(pAng) * &
       sin(lon - lonc)) + cos(lat) * cos(latc) * cos(pAng)
    a23 = -cos(latc) * sin(lat) * cos(lon - lonc) + sin(latc) * cos(lat)
    a31 =  cos(lat) * (sin(latc) * sin(pAng) * cos(lon - lonc) + cos(pAng) * &
       sin(lon - lonc)) - sin(lat) * cos(latc) * sin(pAng)
    a32 = -cos(lat) * (sin(latc) * cos(pAng) * cos(lon - lonc) - sin(pAng) * &
       sin(lon - lonc)) + sin(lat) * cos(latc) * cos(pAng)
    a33 =  cos(lat) * cos(latc) * cos(lon - lonc) + sin(lat) * sin(latc)

    bph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = a11(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) * bxCart(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) +a12(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * byCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) +a13(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * bzCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    bth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = a21(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) * bxCart(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) +a22(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * byCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) +a23(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * bzCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
     br(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = a31(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,ixOmin3:ixOmax3) * bxCart(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,ixOmin3:ixOmax3) +a32(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        ixOmin3:ixOmax3) * byCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        ixOmin3:ixOmax3) +a33(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        ixOmin3:ixOmax3) * bzCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        ixOmin3:ixOmax3)

    A_out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) =         br(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    A_out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2) = -1.0d0*bth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    A_out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3) =        bph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine Cart2SphereVector

  subroutine bipolar_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,A,Bbp)

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! vector potential
    double precision, intent(out)   :: A(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! magnetic field
    double precision, optional, intent(out)   :: Bbp(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)

    double precision :: Aphi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    Aphi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)= q_para*(L_para-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3))/(sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)+d_para)**2)*sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)+d_para)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)-L_para)**2))+q_para*(L_para+x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))/(sqrt(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2+(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)+d_para)**2)*sqrt(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2+(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)+d_para)**2+(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)+L_para)**2))

    A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)=-Aphi(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)/sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)+d_para)**2)
    A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)= 0.d0
    A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)= Aphi(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)+d_para)/sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)+d_para)**2)

    if(present(Bbp)) then
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)+d_para)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)+L_para)**2)**3
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)=       (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)+L_para)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)=                x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)= (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)+d_para)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)+d_para)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)-L_para)**2)**3
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)=      -(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)-L_para)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)=               -x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)=-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)+d_para)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         :)=q_para*Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)
    end if

  end subroutine bipolar_field

  subroutine set_analytic_parker_bipole(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,w)
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: A(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim),Bcart(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim),&
       Bsph(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim),&
       xS(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(:))=0.d0
    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
      call cal_parker_solar_wind(x(ix1,ix2,ix3,1), rc, vs, vout)
      w(ix1,ix2,ix3,rho_)  = (rhob*V_surface)/(Vout*x(ix1,ix2,ix3,1)**2)
      w(ix1,ix2,ix3,mom(1))= vout*w(ix1,ix2,ix3,rho_)
      w(ix1,ix2,ix3,p_)    = w(ix1,ix2,ix3,rho_)*Tiso/(mhd_gamma-1.0d0)
    end do
    end do
    end do

    xS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) * cos(0.50d0*dpi-x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)) * cos(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    xS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) * cos(0.50d0*dpi-x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)) * sin(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    xS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) * sin(0.50d0*dpi-x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
    call bipolar_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xS,A,Bcart)
    call Cart2SphereVector(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,Bcart,Bsph)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mag(:))=Bsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)
  end subroutine set_analytic_parker_bipole

  subroutine boundary_electric_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,fE,s)
    ! specify tangential electric field at physical boundaries 
    ! to fix or drive normal magnetic field
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: qt,qdt
    type(state)                        :: s
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndim:3)

    double precision :: xC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim),xlen1,xlen2,xlen3
    integer :: idir,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixAmin1,&
       ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,ix1,ix2,ix3,ixg1,ixg2,ixg3

    associate(w=>s%w,x=>s%x,ws=>s%ws)

    if(s%is_physical_boundary(1)) then
      ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1;
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      fE(nghostcells,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2:3)=0.0d0
    end if

    ! perfect conductor theta boundaries
    ! is this needed if we want a pole boundary at theta? -------Yihua
    if(block%is_physical_boundary(3)) then
      ixCmin1=ixOmin1-kr(2,1);ixCmin2=ixOmin2-kr(2,2);ixCmin3=ixOmin3-kr(2,3);
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixAmin1=ixCmin1-kr(3,1);ixAmin2=ixCmin2-kr(3,2);ixAmin3=ixCmin3-kr(3,3);
      ixAmax1=ixCmax1;ixAmax2=ixCmax2;ixAmax3=ixCmax3;
      fE(ixAmin1:ixAmax1,ixAmin2,ixAmin3:ixAmax3,1)=0.d0
      ixAmin1=ixCmin1-kr(1,1);ixAmin2=ixCmin2-kr(1,2);ixAmin3=ixCmin3-kr(1,3);
      ixAmax1=ixCmax1;ixAmax2=ixCmax2;ixAmax3=ixCmax3;
      fE(ixAmin1:ixAmax1,ixAmin2,ixAmin3:ixAmax3,3)=0.d0
    end if
    if(block%is_physical_boundary(4)) then
      ixCmin1=ixOmin1-kr(2,1);ixCmin2=ixOmin2-kr(2,2);ixCmin3=ixOmin3-kr(2,3);
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixAmin1=ixCmin1-kr(3,1);ixAmin2=ixCmin2-kr(3,2);ixAmin3=ixCmin3-kr(3,3);
      ixAmax1=ixCmax1;ixAmax2=ixCmax2;ixAmax3=ixCmax3;
      fE(ixAmin1:ixAmax1,ixAmax2,ixAmin3:ixAmax3,1)=0.d0
      ixAmin1=ixCmin1-kr(1,1);ixAmin2=ixCmin2-kr(1,2);ixAmin3=ixCmin3-kr(1,3);
      ixAmax1=ixCmax1;ixAmax2=ixCmax2;ixAmax3=ixCmax3;
      fE(ixAmin1:ixAmax1,ixAmax2,ixAmin3:ixAmax3,3)=0.d0
    end if

    end associate

  end subroutine boundary_electric_field

  subroutine cal_parker_solar_wind(ra_pksw, rc_pksw, vs_pksw, v_pksw)
  use mod_global_parameters
  real*8, intent(in)  :: ra_pksw, rc_pksw, vs_pksw
  real*8, intent(out) :: v_pksw
  real*8              :: W0_pksw, Wn1_pksw, Dr_pksw, nDr_pksw, plw, L_pksw,&
      M_pksw


   Dr_pksw = ((ra_pksw/rc_pksw)**(-4.0d0))*dexp(4.0d0*(1.0d0-&
      (rc_pksw/ra_pksw))-1.0d0)
   nDr_pksw = -1.0d0*Dr_pksw
   if (ra_pksw .le. rc_pksw-0.05d0) then
       if (nDr_pksw .gt. -0.25) then
            W0_pksw = nDr_pksw-nDr_pksw**2.0d0 + &
               1.5d0*nDr_pksw*nDr_pksw*nDr_pksw
       endif
       if (nDr_pksw .le. -0.25) then
           plw=sqrt(2.0d0*(exp(1.0d0)*nDr_pksw+1.0d0))
           W0_pksw=-1.0d0+plw-(1.0d0/3.0d0)*plw**2.0d0+&
              (11.0d0/72.0d0)*plw*plw*plw
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
         M_pksw*(6.0d0-9.0d0*M_pksw+2.0d0*M_pksw*M_pksw)/(&
         6.0d0*L_pksw*L_pksw*L_pksw)
      if (abs(nDr_pksw) .ge. 0.28796090d0) then
         plw=-sqrt(2.0d0*(exp(1.0d0)*nDr_pksw+1.0d0))
         Wn1_pksw=-1.0d0+plw-(1.0d0/3.0d0)*plw*plw+(11.0d0/72.0d0)*plw*plw*plw
      endif
     v_pksw = dsqrt(-1.0d0*vs_pksw*vs_pksw*Wn1_pksw)
   endif
   ! normalize the unit
     v_pksw = v_pksw*1.0d5/unit_velocity  ! Km/s -> normalized velocity

  end subroutine cal_parker_solar_wind

  subroutine specialbound_usr(qt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB,&
        ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    
    double precision :: Qp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: A(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim),Bcart(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim),&
       Bsph(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim),&
       xS(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    integer :: ix1,ix2,ix3,ixOsmin1,ixOsmin2,ixOsmin3,ixOsmax1,ixOsmax2,&
       ixOsmax3,jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3,idir
    double precision :: q,q1,b0,b1,b2
    
    q  = dble(qstretch_baselevel(1))
    q1 = 1.0d0/q
    b0 = -(2.0d0+q1)/(1.0d0+q1)
    b1 = 1.0d0+1.0d0/q1
    b2 = -1.0d0/(1.0d0+q1**2)

    select case(iB)

 case(1)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(:))=0.d0
       do ix3=ixOmin3,ixOmax3
       do ix2=ixOmin2,ixOmax2
       do ix1=ixOmin1,ixOmax1
         call cal_parker_solar_wind(x(ix1,ix2,ix3,1), rc, vs, vout)
         w(ix1,ix2,ix3,rho_)  = (rhob*V_surface)/(Vout*x(ix1,ix2,ix3,1)**2)
         w(ix1,ix2,ix3,mom(1))= vout*w(ix1,ix2,ix3,rho_)
         w(ix1,ix2,ix3,p_)    = w(ix1,ix2,ix3,rho_)*Tiso/(mhd_gamma-1.0d0)
       end do
       end do
       end do

       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
           ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
           ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir)
           ixOsmin3=ixOmin3-kr(3,idir);
           if(stretched_dim(1)) then
             do ix1=ixOsmax1,ixOsmin1,-1
                block%ws(ix1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
                   idir) = ece4(block%level)*block%ws(ix1+4,ixOsmin2:ixOsmax2,&
                   ixOsmin3:ixOsmax3,idir)+ece3(block%level)*block%ws(ix1+3,&
                   ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
                   idir)+ece2(block%level)*block%ws(ix1+2,ixOsmin2:ixOsmax2,&
                   ixOsmin3:ixOsmax3,idir)+ece1(block%level)*block%ws(ix1+1,&
                   ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,idir)
             end do
           else 
             do ix1=ixOsmax1,ixOsmin1,-1
              block%ws(ix1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
                 idir) = 0.12d0*block%ws(ix1+5,ixOsmin2:ixOsmax2,&
                 ixOsmin3:ixOsmax3,idir)-0.76d0*block%ws(ix1+4,&
                 ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
                 idir)+2.08d0*block%ws(ix1+3,ixOsmin2:ixOsmax2,&
                 ixOsmin3:ixOsmax3,idir)-3.36d0*block%ws(ix1+2,&
                 ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
                 idir)+2.92d0*block%ws(ix1+1,ixOsmin2:ixOsmax2,&
                 ixOsmin3:ixOsmax3,idir)

              !block%ws(ix1^%1ixOs^S,idir) = third*&
              !       (-block%ws(ix1+2^%1ixOs^S,idir)&
              !   +4.d0*block%ws(ix1+1^%1ixOs^S,idir))
             end do
            end if 
         end do
         ixOsmin1=ixOmin1-kr(1,1);ixOsmin2=ixOmin2-kr(1,2)
         ixOsmin3=ixOmin3-kr(1,3);ixOsmax1=ixOmax1-kr(1,1)
         ixOsmax2=ixOmax2-kr(1,2);ixOsmax3=ixOmax3-kr(1,3);
         jxOmin1=ixOmin1+nghostcells*kr(1,1)
         jxOmin2=ixOmin2+nghostcells*kr(1,2)
         jxOmin3=ixOmin3+nghostcells*kr(1,3)
         jxOmax1=ixOmax1+nghostcells*kr(1,1)
         jxOmax2=ixOmax2+nghostcells*kr(1,2)
         jxOmax3=ixOmax3+nghostcells*kr(1,3);
         block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
            1)=zero
         do ix1=ixOsmax1,ixOsmin1,-1
           call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Qp)
           block%ws(ix1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,1)=Qp(ix1+1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)*block%dvolume(ix1+1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)/block%surfaceC(ix1,&
              ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,1)
         end do
         call mhd_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
            ixOmax3,block)

         ! Keep the cell-centered ghost magnetic field consistent with the
         ! analytical buried-charge bipole used in the interior initialization.
         xS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            1) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            1) * cos(0.50d0*dpi-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2)) * cos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,3))
         xS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            2) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            1) * cos(0.50d0*dpi-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2)) * sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,3))
         xS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            3) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            1) * sin(0.50d0*dpi-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))
         call bipolar_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,xS,A,Bcart)
         call Cart2SphereVector(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
            ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,Bcart,&
            Bsph)
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(:))=Bsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)
       else
         do ix1=ixOmax1,ixOmin1,-1
           w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))=third* (-w(ix1+2,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)) +4.0d0*w(ix1+1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)))
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
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) = (-b1/b0)*w(ix1-1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)+(-b2/b0)*w(ix1-2,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)   = (-b1/b0)*w(ix1-1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)+(-b2/b0)*w(ix1-2,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1)) = (-b1/b0)*w(ix1-1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))+(-b2/b0)*w(ix1-2,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2:3)) = zero
         end do

       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
             ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
             ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir)
             ixOsmin3=ixOmin3-kr(3,idir);
             do ix1=ixOsmin1,ixOsmax1
                block%ws(ix1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,idir) = 0.0d0
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
         ixOsmin1=ixOmin1;ixOsmin2=ixOmin2;ixOsmin3=ixOmin3;ixOsmax1=ixOmax1
         ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
         jxOmin1=ixOmin1-nghostcells*kr(1,1)
         jxOmin2=ixOmin2-nghostcells*kr(1,2)
         jxOmin3=ixOmin3-nghostcells*kr(1,3)
         jxOmax1=ixOmax1-nghostcells*kr(1,1)
         jxOmax2=ixOmax2-nghostcells*kr(1,2)
         jxOmax3=ixOmax3-nghostcells*kr(1,3);
         block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
            1)=zero
         do ix1=ixOsmin1,ixOsmax1
           call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Qp)
           block%ws(ix1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,1)=-Qp(ix1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)*block%dvolume(ix1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)/block%surfaceC(ix1,&
              ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,1)
         end do
         call mhd_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
            ixOmax3,block)
       else
         do ix1=ixOmin1,ixOmax1
           w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))=third* (-w(ix1-2,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1)) +4.0d0*w(ix1-1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1)))
         enddo
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2:3))=0.d0
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
       call set_analytic_parker_bipole(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,w)

       if(stagger_grid) then
         ! The centered analytical overwrite is fine for ghost cells, but the
         ! CT theta-cut faces still need the divergence-consistent reconstruction.
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
           ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir)
           ixOsmin3=ixOmin3-kr(3,idir);
           do ix2=ixOsmax2,ixOsmin2,-1
             block%ws(ixOsmin1:ixOsmax1,ix2,ixOsmin3:ixOsmax3,&
                idir)=third*(-block%ws(ixOsmin1:ixOsmax1,ix2+2,&
                ixOsmin3:ixOsmax3,idir)+4.d0*block%ws(ixOsmin1:ixOsmax1,ix2+1,&
                ixOsmin3:ixOsmax3,idir))
           end do
         end do
         ixOsmin1=ixOmin1-kr(2,1);ixOsmin2=ixOmin2-kr(2,2)
         ixOsmin3=ixOmin3-kr(2,3);ixOsmax1=ixOmax1-kr(2,1)
         ixOsmax2=ixOmax2-kr(2,2);ixOsmax3=ixOmax3-kr(2,3);
         jxOmin1=ixOmin1+nghostcells*kr(2,1)
         jxOmin2=ixOmin2+nghostcells*kr(2,2)
         jxOmin3=ixOmin3+nghostcells*kr(2,3)
         jxOmax1=ixOmax1+nghostcells*kr(2,1)
         jxOmax2=ixOmax2+nghostcells*kr(2,2)
         jxOmax3=ixOmax3+nghostcells*kr(2,3);
         block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
            2)=zero
         do ix2=ixOsmax2,ixOsmin2,-1
           call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Qp)
           block%ws(ixOsmin1:ixOsmax1,ix2,ixOsmin3:ixOsmax3,&
              2)=Qp(ixOmin1:ixOmax1,ix2+1,&
              ixOmin3:ixOmax3)*block%dvolume(ixOmin1:ixOmax1,ix2+1,&
              ixOmin3:ixOmax3)/block%surfaceC(ixOsmin1:ixOsmax1,ix2,&
              ixOsmin3:ixOsmax3,2)
         end do
         call mhd_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
            ixOmax3,block)
         call set_analytic_parker_bipole(ixImin1,ixImin2,ixImin3,ixImax1,&
            ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,&
            w)
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
       call set_analytic_parker_bipole(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,w)

       if(stagger_grid) then
         ! The centered analytical overwrite is fine for ghost cells, but the
         ! CT theta-cut faces still need the divergence-consistent reconstruction.
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
           ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir)
           ixOsmin3=ixOmin3-kr(3,idir);
           do ix2=ixOsmin2,ixOsmax2
             block%ws(ixOsmin1:ixOsmax1,ix2,ixOsmin3:ixOsmax3,&
                idir)=third*(-block%ws(ixOsmin1:ixOsmax1,ix2-2,&
                ixOsmin3:ixOsmax3,idir)+4.d0*block%ws(ixOsmin1:ixOsmax1,ix2-1,&
                ixOsmin3:ixOsmax3,idir))
           end do
         end do
         ixOsmin1=ixOmin1;ixOsmin2=ixOmin2;ixOsmin3=ixOmin3;ixOsmax1=ixOmax1
         ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
         jxOmin1=ixOmin1-nghostcells*kr(2,1)
         jxOmin2=ixOmin2-nghostcells*kr(2,2)
         jxOmin3=ixOmin3-nghostcells*kr(2,3)
         jxOmax1=ixOmax1-nghostcells*kr(2,1)
         jxOmax2=ixOmax2-nghostcells*kr(2,2)
         jxOmax3=ixOmax3-nghostcells*kr(2,3);
         block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
            2)=zero
         do ix2=ixOsmin2,ixOsmax2
           call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Qp)
           block%ws(ixOsmin1:ixOsmax1,ix2,ixOsmin3:ixOsmax3,&
              2)=-Qp(ixOmin1:ixOmax1,ix2,&
              ixOmin3:ixOmax3)*block%dvolume(ixOmin1:ixOmax1,ix2,&
              ixOmin3:ixOmax3)/block%surfaceC(ixOsmin1:ixOsmax1,ix2,&
              ixOsmin3:ixOsmax3,2)
         end do
         call mhd_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
            ixOmax3,block)
         call set_analytic_parker_bipole(ixImin1,ixImin2,ixImin3,ixImax1,&
            ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,&
            w)
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
  subroutine getggrav(ggrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: ggrid(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    ggrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=usr_grav*SRadius**2/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1))**2
  end subroutine

  !==============================================================================
  ! Purpose: get gravity field
  !==============================================================================
  subroutine gravity(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,ndim)
    double precision                :: ggrid(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    gravity_field=0.d0
    call getggrav(ggrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x)
    gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)=ggrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
  end subroutine gravity

  !==============================================================================
  ! Purpose: Enforce additional refinement or coarsening. One can use the
  !          coordinate info in x and/or time qt=t_n and w(t_n) values w.
  !==============================================================================
  subroutine special_refine_grid(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x,&
     refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qt, w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision :: rmin, thetamin, thetamax
    double precision :: dr_block, dtheta_block, dphi_block
    double precision :: jb_eta_max, jb_eta_mean
    double precision :: bmag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        current_mag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        delta_min(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        jb_eta(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,3)
    logical :: safe_theta, force_inner, force_sheet
    integer :: idirmin

    refine=0
    coarsen=0

    rmin=minval(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))
    thetamin=minval(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
    thetamax=maxval(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))

    safe_theta = (thetamin >= xprobmin2 + amr_theta_guard) .and. (thetamax <= &
       xprobmax2 - amr_theta_guard)

    dr_block=(maxval(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1))-minval(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)))/dble(max(1,ixOmax1-ixOmin1))
    dtheta_block=(maxval(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2))-minval(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)))/dble(max(1,ixOmax2-ixOmin2))
    dphi_block=(maxval(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3))-minval(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)))/dble(max(1,ixOmax3-ixOmin3))

    dxlevel(1)=rnode(rpdx1_,igrid)
    dxlevel(2)=rnode(rpdx2_,igrid)
    dxlevel(3)=rnode(rpdx3_,igrid)
    block=>ps(igrid)

    bmag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)=sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       mag(1))**2 + w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       mag(2))**2 + w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       mag(3))**2)
    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    current_mag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)=sqrt(current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1)**2 + current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,2)**2 + current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,3)**2)
    delta_min(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=min(dr_block,&
        x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1)*dtheta_block,&
        x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1)*max(dsin(x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2)),&
       smalldouble)*dphi_block)
    jb_eta(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)=current_mag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)*delta_min(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)/(bmag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)+smalldouble)
    jb_eta_max=maxval(jb_eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    jb_eta_mean=sum(jb_eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))/dble(size(jb_eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)))

    force_inner = safe_theta .and. (rmin <= amr_r_core_max)
    force_sheet = safe_theta .and. (rmin <= amr_r_sheet_max) .and. (jb_eta_max &
       >= amr_jb_eta_max_min) .and. (jb_eta_mean >= amr_jb_eta_mean_min)

    if(level < amr_force_max_level .and. (force_inner .or. force_sheet)) then
      refine=1
      coarsen=-1
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
  subroutine specialvar_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision                   :: ws1(ixGslo1:ixGshi1,ixGslo2:ixGshi2,&
       ixGslo3:ixGshi3,1:ndim)
    double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw+nwauxio)
    double precision                   :: w_om(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision                   :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision                   :: bvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision                   :: patchwi(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision                   :: xOmin1,xOmin2,xOmin3,xOmax1,xOmax2,&
       xOmax3
    double precision                   :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: csound2(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)

    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,3)
    double precision :: xlen1,xlen2,xlen3
    integer :: ix1,ix2,ix3,ixbc1,ixbc2,ixbc3,idirmin,idir,jdir,kdir,hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,idim

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+2)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+3)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))
    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+4)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+5)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+6)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
    xOmin1 = xprobmin1 + 0.05d0*(xprobmax1-xprobmin1)
    xOmin2 = xprobmin2 + 0.05d0*(xprobmax2-xprobmin2)
    xOmin3 = xprobmin3 + 0.05d0*(xprobmax3-xprobmin3)
    xOmax1 = xprobmax1 - 0.05d0*(xprobmax1-xprobmin1)
    xOmax2 = xprobmax2 - 0.05d0*(xprobmax2-xprobmin2)
    xOmax3 = xprobmax3 - 0.05d0*(xprobmax3-xprobmin3)
    xOmin1 = xprobmin1

    qvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)=zero
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)=zero
    tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=zero
    tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=zero
    w_om=zero
    patchwi=zero

    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
        if( x(ix1,ix2,ix3,1) > xOmin1 .and. x(ix1,ix2,ix3,&
           1) < xOmax1  .and.  x(ix1,ix2,ix3,2) > xOmin2 .and. x(ix1,ix2,ix3,&
           2) < xOmax2  .and.  x(ix1,ix2,ix3,3) > xOmin3 .and. x(ix1,ix2,ix3,&
           3) < xOmax3 ) then
          patchwi(ix1,ix2,ix3)=1.0d0
        else
          patchwi(ix1,ix2,ix3)=0.0d0
        endif
    end do
    end do
    end do

    if(B0field) then
      bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         :)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(:))+block%b0(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(:),0)
    else
      bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         :)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(:))
    endif
    do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
       if(lvc(idir,jdir,kdir)/=0)then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,jdir)*bvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,kdir)
          if(lvc(idir,jdir,kdir)==1)then
             qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                idir)=qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                idir)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          else
             qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                idir)=qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                idir)-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          endif
       endif
    enddo; enddo; enddo

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+7)=patchwi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*sqrt(sum(qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,:)**2,ndim+1))/sqrt(sum(bvec(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)**2,ndim+1))
    w_om(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:nw)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nw)
    call get_normalized_divb(w_om,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp1)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+8)=tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+9)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))*unit_velocity/(w(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)*1.0d5)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+10)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))*unit_velocity/(w(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)*1.0d5)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+11)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))*unit_velocity/(w(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)*1.0d5)

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
    
    double precision :: theta_zone, outer_zone, tol_add, dist_theta_cut

    tol_add=0.d0
    theta_zone=max(amr_theta_guard, amr_theta_relax_frac*(xprobmax2-&
       xprobmin2))
    outer_zone=amr_outer_relax_frac*(xprobmax1-xprobmin1)
    dist_theta_cut=min(xlocal(2)-xprobmin2, xprobmax2-xlocal(2))

    if(dist_theta_cut < theta_zone) then
      tol_add=tol_add + (1.d0-dist_theta_cut/theta_zone)*0.5d0
    endif

    if(xprobmax1-xlocal(1) < outer_zone) then
      tol_add=tol_add + (1.d0-(xprobmax1-xlocal(1))/outer_zone)*0.5d0
    endif

    tolerance=tolerance+tol_add

  end subroutine specialthreshold

end module mod_usr
