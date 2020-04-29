!> Functionality to support the different types of grid and abstraction between global grids
!! and local ones after decomposition
!!
!! Currently MONC supports the Arakawa C grid
module grids_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif

  !> Grid index parameters
  integer, parameter, public :: Z_INDEX = 1, Y_INDEX = 2, X_INDEX = 3
  !> Grid type parameters (usually applied to each dimension of a prognostic)
  integer, parameter, public :: PRIMAL_GRID=1, DUAL_GRID=2

  !> The configuration of the grid horizontally
  type, public :: horizontal_grid_configuration_type
     real(kind=DEFAULT_PRECISION) dx, dy     !< Grid spacings in x and y directions (m)
     real(kind=DEFAULT_PRECISION) cx, cy     !< Reciprocal of the horizontal grid spacings
     real(kind=DEFAULT_PRECISION) cx2, cy2   !< Reciprocal of the square of the grid spacings
     real(kind=DEFAULT_PRECISION) cxy        !< 1/(DX*DY)
     real(kind=DEFAULT_PRECISION) tcx, tcy
  end type horizontal_grid_configuration_type

  !> The configuration of the grid vertically
  type, public :: vertical_grid_configuration_type
     real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: &
          z,&      !< Heights at w levels (m)
          zn,&     !< Heights at pressure levels (m)
          dz,&     !< Vertical spacing between K and K-1 w levels
          dzn,&    !< Vertical spacing between K and K-1 p levels
          rdz,&    !< Reciprocal of DZ
          rdzn,&   !< Reciprocal of DZN
          rho,&    !< Density at w levels (kg/m3)
          rhon,&   !< Density at p levels (kg/m3)
          thref,&  !< Reference potential temperature (K)
          dthref,& !< Gradient of thref (K)
          tref,&   !< Reference temperature (K)
          theta_init,& !< Initial profile of potential temperature
          temp_init,&  !< Initial profile of absolute temperature
          rh_init, &   !< Initial profile of relative humidity
          u_init,&     !< Initial profile of u
          v_init,&     !< Initial profile of v
          prefn,& !< Reference pressure (Pa)
          pdiff,& !< Difference between pressure levels (Pa)
          prefrcp,&
          rprefrcp,&
          czb,&   !< (rho(k-1)/rhon(k))/(dz(k)*dzn(k)) use for diffusion onto p level from below
          cza,&   !< (rho(k)/rhon(k))/(dz(k)*dzn(k+1)) use for diffusion onto p level from above
          czg,&   !< CZB-CZA for tridiagonal solver in POISSON
          czh,&   !< CZB*CZA for tridiagonal solver in POISSON      
          tzc1,&  !< 0.25*rdz(k)*rho(k-1)/rhon(k) for advection onto p-level from below
          tzc2,&  !< 0.25*rdz(k)*rho(k)/rhon(k) for advection onto p-level from above
          tzd1,&  !< 0.25*rdzn(k+1)*rhon(k)/rho(k) for advection onto w-level from below
          tzd2,&  !< 0.25*rdzn(k+1)*rhon(k+1)/rho(k) for advection onto w-level from above 
          w_subs,& !< Subsidence velocity
          olubar,& !< Current U mean
          savolubar,&          
          olvbar,& !< Current V mean
          savolvbar,&
          olthbar,& !< Current theta mean
          olzubar,& !< Previous timestep U mean
          olzvbar,& !< Previous timestep V mean
          olzthbar,& !< Previous theta mean
          dmpco,& !< Damping coefficient on pressure levels
          dmpcoz,& !< Damping coefficient on w-levels
          tstarpr,& ! Temperature about which Taylor Expansion
          qsat,&
          dqsatdt,&
          qsatfac,&
          rneutml,&
          rneutml_sq,&
          buoy_co, &
          theta_rand, &!< profile of amplitude of theta perturbations
          theta_force, & !<profile of forcing term for theta
          u_force, & !<profile of forcing term for u
          v_force, & !<profile of forcing term for v
          w_up,    & !<profile of updraft threshold
          w_dwn,   & !<profile of downdraft threshold
          w_rand     !<profile of amplitude of w perturbations



     real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_init !< Initial profile of q variables
     real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_rand !< Initial profile of amplitude of q perturbations
     real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_force !< Profiles of forcing terms for q variables

     real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: olqbar,olzqbar
    ! time varying forcing terms
     real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: wsubs_time_vary 
     
     real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: hgd
     real(kind=DEFAULT_PRECISION) :: czn, zlogm, zlogth, vk_on_zlogm
     integer, dimension(:), allocatable :: kgd
     integer :: kdmpmin
  end type vertical_grid_configuration_type

  !> Wraps the dimensional configuration types
  type, public :: grid_configuration_type
     type(horizontal_grid_configuration_type) :: horizontal !< Configuration horizontally
     type(vertical_grid_configuration_type) :: vertical !< Configuration vertically
  end type grid_configuration_type

  !> Defines the global grid
  type, public :: global_grid_type
     type(grid_configuration_type) :: configuration !< Configuration of the grid
     real(kind=DEFAULT_PRECISION), dimension(3) :: bottom,&     !< Bottom (lowest) bounds in each dimension
          top,&        !< Top (highest) bounds in each dimension
          resolution   !< The resolution of the grid in each dimension
     integer, dimension(3) ::    size         !< Number of grid points in each dimension
     logical, dimension(3) :: active = (/ .false., .false., .false. /) !< Whether a specific dimension is active
     integer :: dimensions = 0 !< Number of active dimensions
  end type global_grid_type

  !> Defined the local grid, i.e. the grid held on this process after decomposition
  type, public :: local_grid_type
     integer, dimension(3) :: start,& !< Global start coordinates of local grid
          end,& !< Global end coordinates of local grid
          size,& !< Grid points held locally
          halo_size,& !< Grid point size of the halo (halo_depth)
          local_domain_start_index,& !< The start index of computation data (local data is halo->data->halo so this precomputes the data start)
          local_domain_end_index !< The end index of computation data (local data is halo->data->halo so this precomputes the data end)
     logical, dimension(3) :: active = (/ .false., .false., .false. /) !< Dimensions which are active
     integer, dimension(:,:), allocatable :: neighbours , corner_neighbours !< Neighbouring process Id per dimension
     integer :: dimensions = 0 !< Number of active dimensions
  end type local_grid_type
end module grids_mod
