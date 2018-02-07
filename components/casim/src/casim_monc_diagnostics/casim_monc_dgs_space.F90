module casim_monc_dgs_space

  use generic_diagnostic_variables, only: diaglist
  use state_mod, only :  model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION
  use mphys_switches, only: l_warm

  implicit none

  ! Method:
  ! Allocates and de-allocates the diagnostics space for
  ! casim monc diags. This is different to the allocation in
  ! in casdiags struct (generic_diagnostics) since the
  ! this declares the size of the local decomposition
  ! which needs to be populated for the IO server

  TYPE casim_monc_dglist

     !--------------------------------
     ! 2D variable arrays
     !--------------------------------
     ! Surface Precipitation rates
     REAL, ALLOCATABLE :: precip(:,:)
     REAL, ALLOCATABLE :: SurfaceRainR(:,:)
     REAL, ALLOCATABLE :: SurfaceSnowR(:,:)
     REAL, ALLOCATABLE :: SurfaceGraupR(:,:)

     ! Process rate diagnostics
     REAL, ALLOCATABLE :: phomc(:,:,:)
     REAL, ALLOCATABLE :: pinuc(:,:,:)
     REAL, ALLOCATABLE :: pidep(:,:,:)
     REAL, ALLOCATABLE :: psdep(:,:,:)
     REAL, ALLOCATABLE :: piacw(:,:,:)
     REAL, ALLOCATABLE :: psacw(:,:,:)
     REAL, ALLOCATABLE :: psacr(:,:,:)
     REAL, ALLOCATABLE :: pisub(:,:,:)
     REAL, ALLOCATABLE :: pssub(:,:,:)
     REAL, ALLOCATABLE :: pimlt(:,:,:)
     REAL, ALLOCATABLE :: psmlt(:,:,:)
     REAL, ALLOCATABLE :: psaut(:,:,:)
     REAL, ALLOCATABLE :: psaci(:,:,:)
     REAL, ALLOCATABLE :: praut(:,:,:)
     REAL, ALLOCATABLE :: pracw(:,:,:)
     REAL, ALLOCATABLE :: prevp(:,:,:)
     REAL, ALLOCATABLE :: pgacw(:,:,:)
     REAL, ALLOCATABLE :: pgacs(:,:,:)
     REAL, ALLOCATABLE :: pgmlt(:,:,:)
     REAL, ALLOCATABLE :: pgsub(:,:,:)
     REAL, ALLOCATABLE :: psedi(:,:,:)
     REAL, ALLOCATABLE :: pseds(:,:,:)
     REAL, ALLOCATABLE :: psedr(:,:,:)
     REAL, ALLOCATABLE :: psedg(:,:,:)
     REAL, ALLOCATABLE :: psedl(:,:,:)
     REAL, ALLOCATABLE :: pcond(:,:,:)

     !---------------------------------
     ! 3D variables logical for
     ! theta tendencies (based on LEM and
     ! MONC, should work with UM)
     !--------------------------------
     REAL, ALLOCATABLE :: dth_total(:,:,:)
     REAL, ALLOCATABLE :: dth_cond_evap(:,:,:)
     
     !---------------------------------
     ! 3D variables for
     ! mass tendencies (based on LEM and
     ! MONC, should work with UM)
     !--------------------------------
     REAL, ALLOCATABLE :: dqv_total(:,:,:)
     REAL, ALLOCATABLE :: dqv_cond_evap(:,:,:)
     REAL, ALLOCATABLE :: dqc(:,:,:)
     REAL, ALLOCATABLE :: dqr(:,:,:)
     REAL, ALLOCATABLE :: dqi(:,:,:)
     REAL, ALLOCATABLE :: dqs(:,:,:)
     REAL, ALLOCATABLE :: dqg(:,:,:)
     
  end type casim_monc_dglist

  type (casim_monc_dglist) :: casim_monc_dgs

contains

  subroutine allocate_casim_monc_dgs_space(current_state, casdiags)
    
    type(model_state_type), target, intent(inout) :: current_state
    type(diaglist), target, intent(in) :: casdiags
    
    integer :: z_size, y_size_local, x_size_local

    z_size = current_state%local_grid%size(Z_INDEX)
    y_size_local = current_state%local_grid%size(Y_INDEX)
    x_size_local = current_state%local_grid%size(X_INDEX)

    ! precipitation arrays
    allocate ( casim_monc_dgs % precip(y_size_local, x_size_local) )
    casim_monc_dgs % precip(:,:) = 0.0_DEFAULT_PRECISION

    allocate ( casim_monc_dgs % SurfaceRainR(y_size_local, x_size_local) )
    casim_monc_dgs % SurfaceRainR(:, :) = 0.0_DEFAULT_PRECISION
    
    ! Process rate diagnostics
       
    allocate ( casim_monc_dgs % psedl(z_size, y_size_local, x_size_local) )
    casim_monc_dgs % psedl(:, :, :) = 0.0_DEFAULT_PRECISION
    
    allocate ( casim_monc_dgs % pcond(z_size, y_size_local, x_size_local) )
    casim_monc_dgs % pcond(:, :, :) = 0.0_DEFAULT_PRECISION
       
    allocate ( casim_monc_dgs % praut(z_size, y_size_local, x_size_local) )
    casim_monc_dgs % praut(:, :, :) = 0.0_DEFAULT_PRECISION
    
    allocate ( casim_monc_dgs % pracw(z_size, y_size_local, x_size_local) )
    casim_monc_dgs % pracw(:, :, :) = 0.0_DEFAULT_PRECISION
    
    allocate ( casim_monc_dgs % prevp(z_size, y_size_local, x_size_local) )
    casim_monc_dgs % prevp(:, :, :) = 0.0_DEFAULT_PRECISION
    
    allocate ( casim_monc_dgs % psedr(z_size, y_size_local, x_size_local) )
    casim_monc_dgs % psedr(:, :, :) = 0.0_DEFAULT_PRECISION

    ! potential temperature and mass tendencies
    allocate ( casim_monc_dgs % dth_total(z_size, y_size_local, x_size_local) )
    casim_monc_dgs % dth_total(:, :, :) = 0.0_DEFAULT_PRECISION
    
    allocate ( casim_monc_dgs % dth_cond_evap(z_size, y_size_local, x_size_local) )
    casim_monc_dgs % dth_cond_evap(:, :, :) = 0.0_DEFAULT_PRECISION

    allocate ( casim_monc_dgs % dqv_total(z_size, y_size_local, x_size_local) )
    casim_monc_dgs % dqv_total(:, :, :) = 0.0_DEFAULT_PRECISION

    allocate ( casim_monc_dgs % dqv_cond_evap(z_size, y_size_local, x_size_local) )
    casim_monc_dgs % dqv_cond_evap(:, :, :) = 0.0_DEFAULT_PRECISION

    allocate ( casim_monc_dgs % dqc(z_size, y_size_local, x_size_local) )
    casim_monc_dgs % dqc(:, :, :) = 0.0_DEFAULT_PRECISION
    
    allocate ( casim_monc_dgs % dqr(z_size, y_size_local, x_size_local) )
    casim_monc_dgs % dqr(:, :, :) = 0.0_DEFAULT_PRECISION
    
    if (.not. l_warm) then

       ! ice/cold process rates
       allocate ( casim_monc_dgs % SurfaceSnowR(y_size_local, x_size_local) )
       casim_monc_dgs % SurfaceSnowR(:,:) =  0.0_DEFAULT_PRECISION

       allocate ( casim_monc_dgs % SurfaceGraupR(y_size_local, x_size_local) )
       casim_monc_dgs % SurfaceGraupR(:,:) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % phomc(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % phomc(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % pinuc(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % pinuc(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % pidep(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % pidep(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % piacw(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % piacw(:, :, :) = 0.0_DEFAULT_PRECISION

       allocate ( casim_monc_dgs % pisub(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % pisub(:, :, :) = 0.0_DEFAULT_PRECISION

       allocate ( casim_monc_dgs % pimlt(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % pimlt(:, :, :) = 0.0_DEFAULT_PRECISION

       allocate ( casim_monc_dgs % psedi(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % psedi(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % psdep(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % psdep(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % psacw(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % psacw(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % psacr(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % psacr(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % pssub(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % pssub(:, :, :) = 0.0_DEFAULT_PRECISION
    
       allocate ( casim_monc_dgs % psmlt(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % psmlt(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % psaut(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % psaut(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % psaci(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % psaci(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % pgacw(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % pgacw(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % pgacs(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % pgacs(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % pgmlt(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % pgmlt(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % pgsub(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % pgsub(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % pseds(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % pseds(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % psedg(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % psedg(:, :, :) = 0.0_DEFAULT_PRECISION

       ! ice, snow, graupel tendencies
       allocate ( casim_monc_dgs % dqi(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % dqi(:, :, :) = 0.0_DEFAULT_PRECISION
       
       allocate ( casim_monc_dgs % dqs(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % dqs(:, :, :) = 0.0_DEFAULT_PRECISION

       allocate ( casim_monc_dgs % dqg(z_size, y_size_local, x_size_local) )
       casim_monc_dgs % dqg(:, :, :) = 0.0_DEFAULT_PRECISION
       
    endif
    
  end subroutine allocate_casim_monc_dgs_space

  subroutine populate_casim_monc_dg(current_state, casdiags )
    
    type(model_state_type), target, intent(in) :: current_state
    type(diaglist), target, intent(in) :: casdiags

    INTEGER :: icol, jcol, target_x_index, target_y_index, k

    icol=current_state%column_local_x
    jcol=current_state%column_local_y
    target_y_index=jcol-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=icol-current_state%local_grid%halo_size(X_INDEX)

    casim_monc_dgs % precip(target_y_index,target_x_index) = &
         casdiags % SurfaceRainR(1,1)
    
    casim_monc_dgs % SurfaceRainR(target_y_index,target_x_index) = &
         casdiags % SurfaceRainR(1,1)
    
    casim_monc_dgs % pcond(:,target_y_index,target_x_index) = &
         casdiags % pcond(1,1,:)
    casim_monc_dgs % psedl(:,target_y_index,target_x_index) = &
         casdiags % psedl(1,1,:)
    casim_monc_dgs % praut(:,target_y_index,target_x_index) = &
         casdiags % praut(1,1,:)
    casim_monc_dgs % pracw(:,target_y_index,target_x_index) = &
         casdiags % praut(1,1,:)
    casim_monc_dgs % prevp(:,target_y_index,target_x_index) = &
         casdiags % prevp(1,1,:)
    casim_monc_dgs % psedr(:,target_y_index,target_x_index) = &
         casdiags % psedr(1,1,:)

    ! potential temperature and mass tendencies
    casim_monc_dgs % dth_total(:,target_y_index,target_x_index) = &
         casdiags % dth_total(1,1,:)
    casim_monc_dgs % dth_cond_evap(:,target_y_index,target_x_index) = &
         casdiags % dth_cond_evap(1,1,:)
    casim_monc_dgs % dqv_total(:,target_y_index,target_x_index) = &
         casdiags % dqv_total(1,1,:)
    casim_monc_dgs % dqv_cond_evap(:,target_y_index,target_x_index) = &
         casdiags % dqv_cond_evap(1,1,:)
    casim_monc_dgs % dqc(:,target_y_index,target_x_index) = &
         casdiags % dqc(1,1,:)
    casim_monc_dgs % dqr(:,target_y_index,target_x_index) = &
         casdiags % dqr(1,1,:)
    
    if (.not. l_warm) then
       casim_monc_dgs % precip(target_y_index,target_x_index) = &
            casdiags % SurfaceRainR(1,1) + casdiags % SurfaceSnowR(1,1)
       casim_monc_dgs % SurfaceSnowR(target_y_index,target_x_index) = &
            casdiags % SurfaceSnowR(1,1)
       casim_monc_dgs % SurfaceGraupR(target_y_index,target_x_index) = &
            casdiags % SurfaceGraupR(1,1)
       casim_monc_dgs % phomc(:,target_y_index,target_x_index) = &
            casdiags % psedr(1,1,:) 
       casim_monc_dgs % pinuc(:,target_y_index,target_x_index) = &
            casdiags % pinuc(1,1,:)
       casim_monc_dgs % pidep(:,target_y_index,target_x_index) = &
            casdiags % pidep(1,1,:)
       casim_monc_dgs % piacw(:,target_y_index,target_x_index) = &
            casdiags % piacw(1,1,:)
       casim_monc_dgs % pisub(:,target_y_index,target_x_index) = &
            casdiags % pisub(1,1,:)
       casim_monc_dgs % pimlt(:,target_y_index,target_x_index) = &
            casdiags % pimlt(1,1,:)
       casim_monc_dgs % psedi(:,target_y_index,target_x_index) = &
            casdiags % psedi(1,1,:)  
       casim_monc_dgs % psmlt(:,target_y_index,target_x_index) = &
            casdiags % psmlt(1,1,:)
       casim_monc_dgs % psaut(:,target_y_index,target_x_index) = &
            casdiags % psaut(1,1,:) 
       casim_monc_dgs % psaci(:,target_y_index,target_x_index) = &
            casdiags % psaci(1,1,:)
       casim_monc_dgs % psacw(:,target_y_index,target_x_index) = &
            casdiags % psacw(1,1,:)
       casim_monc_dgs % psacr(:,target_y_index,target_x_index) = &
            casdiags % psacr(1,1,:) 
       casim_monc_dgs % pssub(:,target_y_index,target_x_index) = &
            casdiags % pssub(1,1,:)
       casim_monc_dgs % psdep(:,target_y_index,target_x_index) = &
            casdiags % psdep(1,1,:)
       casim_monc_dgs % pseds(:,target_y_index,target_x_index) = &
            casdiags % pseds(1,1,:)
       casim_monc_dgs % pgacw(:,target_y_index,target_x_index) = &
            casdiags % pgacw(1,1,:)  
       casim_monc_dgs % pgacs(:,target_y_index,target_x_index) = &
            casdiags % pgacs(1,1,:)
       casim_monc_dgs % pgmlt(:,target_y_index,target_x_index) = &
            casdiags % pgmlt(1,1,:)
       casim_monc_dgs % pgsub(:,target_y_index,target_x_index) = &
            casdiags % pgsub(1,1,:)
       casim_monc_dgs % psedg(:,target_y_index,target_x_index) = &
            casdiags % psedg(1,1,:)

       casim_monc_dgs % dqi(:,target_y_index,target_x_index) = &
            casdiags % dqi(1,1,:)
       casim_monc_dgs % dqs(:,target_y_index,target_x_index) = &
            casdiags % dqs(1,1,:)
       casim_monc_dgs % dqg(:,target_y_index,target_x_index) = &
            casdiags % dqg(1,1,:)
       
    endif
    
  end subroutine populate_casim_monc_dg
  
end module casim_monc_dgs_space

