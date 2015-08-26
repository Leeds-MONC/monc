module diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: total_points, hqlmax, ncl_col
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tempfac
  real(kind=DEFAULT_PRECISION) :: totqv, totql, totqz, wmax, wmin, qlmax, qlcrit, cltop_av, clbas_av, cltop, clbas

  public diagnostics_get_descriptor

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function diagnostics_get_descriptor()
    diagnostics_get_descriptor%name="diagnostics"
    diagnostics_get_descriptor%version=0.1

    diagnostics_get_descriptor%initialisation=>initialisation_callback
    diagnostics_get_descriptor%timestep=>timestep_callback

    diagnostics_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    diagnostics_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(diagnostics_get_descriptor%published_fields(14))

    diagnostics_get_descriptor%published_fields(1)="totpnts"
    diagnostics_get_descriptor%published_fields(2)="zn_top"
    diagnostics_get_descriptor%published_fields(3)="wmax_local"
    diagnostics_get_descriptor%published_fields(4)="wmin_local"
    diagnostics_get_descriptor%published_fields(5)="totqv_local"
    diagnostics_get_descriptor%published_fields(6)="totql_local"
    diagnostics_get_descriptor%published_fields(7)="totqz_local"
    diagnostics_get_descriptor%published_fields(8)="qlmax_local"
    diagnostics_get_descriptor%published_fields(9)="hqlmax_local"
    diagnostics_get_descriptor%published_fields(10)="cltop_local"
    diagnostics_get_descriptor%published_fields(11)="clbas_local"
    diagnostics_get_descriptor%published_fields(12)="cltop_av_local"
    diagnostics_get_descriptor%published_fields(13)="clbas_av_local"
    diagnostics_get_descriptor%published_fields(14)="ncl_col_local"
  end function diagnostics_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    qlcrit=options_get_real(current_state%options_database, "qlcrit")

    total_points=current_state%global_grid%size(Z_INDEX) * current_state%global_grid%size(Y_INDEX) * &
           current_state%global_grid%size(X_INDEX)

    allocate(tempfac(current_state%local_grid%size(Z_INDEX)))
    do k=2, current_state%local_grid%size(Z_INDEX)
      tempfac(k)=current_state%global_grid%configuration%vertical%dz(k)*&
           current_state%global_grid%configuration%vertical%rhon(k)/total_points
    end do    
  end subroutine initialisation_callback  

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION) :: cltop_col, clbas_col

    if (current_state%first_timestep_column) then
      totqv=0.0_DEFAULT_PRECISION
      totql=0.0_DEFAULT_PRECISION
      totqz=0.0_DEFAULT_PRECISION
      wmax=0.0_DEFAULT_PRECISION
      wmin=0.0_DEFAULT_PRECISION
      qlmax=0.0_DEFAULT_PRECISION
      cltop_av=0.0_DEFAULT_PRECISION
      clbas_av=0.0_DEFAULT_PRECISION
      cltop=0.0_DEFAULT_PRECISION
      clbas=0.0_DEFAULT_PRECISION
      ncl_col=0
    end if
    if (.not. current_state%halo_column) then
      cltop_col=0.0_DEFAULT_PRECISION
      clbas_col=0.0_DEFAULT_PRECISION
      do k=2, current_state%local_grid%size(Z_INDEX)
        wmax=max(wmax, current_state%w%data(k, current_state%column_local_y, current_state%column_local_x))
        wmin=min(wmin, current_state%w%data(k, current_state%column_local_y, current_state%column_local_x))
        if (qlmax .lt. current_state%q(current_state%liquid_water_mixing_ratio_index)%data(k, &
             current_state%column_local_y, current_state%column_local_x)) then
          qlmax=max(qlmax, current_state%q(current_state%liquid_water_mixing_ratio_index)%data(k, &
               current_state%column_local_y, current_state%column_local_x))
          hqlmax=k
        end if

        if (current_state%q(current_state%liquid_water_mixing_ratio_index)%data(k, &
             current_state%column_local_y, current_state%column_local_x) .gt. qlcrit) then
          cltop_col=current_state%global_grid%configuration%vertical%zn(k)
          cltop=max(cltop, current_state%global_grid%configuration%vertical%zn(k))
          clbas=min(clbas, current_state%global_grid%configuration%vertical%zn(k))
        end if

        if (current_state%q(current_state%liquid_water_mixing_ratio_index)%data(current_state%local_grid%size(Z_INDEX)+1-k, &
             current_state%column_local_y, current_state%column_local_x) .gt. qlcrit) then
          clbas_col=current_state%global_grid%configuration%vertical%zn(current_state%local_grid%size(Z_INDEX)+1-k)
        end if

        totqv=totqv+tempfac(k)*current_state%q(current_state%water_vapour_mixing_ratio_index)%data(k, &
             current_state%column_local_y, current_state%column_local_x)
        totql=totql+tempfac(k)*current_state%q(current_state%liquid_water_mixing_ratio_index)%data(k, &
             current_state%column_local_y, current_state%column_local_x)

        !tql_loc(J,I)=tql_loc(J,I)+DZ(K)*RHON(K)*Q(J,K,IQL,current_state%liquid_water_mixing_ratio_index)
        totqz=totqz+tempfac(k)*current_state%q(current_state%liquid_water_mixing_ratio_index)%data(k, &
             current_state%column_local_y, current_state%column_local_x)*current_state%global_grid%configuration%vertical%zn(k)
      end do
      if (cltop_col .gt. 0.0_DEFAULT_PRECISION) ncl_col=ncl_col+1
      cltop_av=cltop_av+cltop_col
      clbas_av=clbas_av+clbas_col
    end if    
  end subroutine timestep_callback  

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information

    field_information%field_type=COMPONENT_SCALAR_FIELD_TYPE
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    field_information%enabled=.true.

    if (name .eq. "totpnts" .or. name .eq. "hqlmax_local" .or. name .eq. "ncl_col_local") then
      field_information%data_type=COMPONENT_INTEGER_DATA_TYPE
    else
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    end if   
  end subroutine field_information_retrieval_callback

  !> Field value retrieval callback, this returns the value of a specific published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve the value for
  !! @param field_value Populated with the value of the field
  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value    

    if (name .eq. "totpnts") then
      field_value%scalar_int=total_points
    else if (name .eq. "zn_top") then
      field_value%scalar_real=current_state%global_grid%configuration%vertical%zn(current_state%global_grid%size(Z_INDEX))
    else if (name .eq. "wmax_local") then
      field_value%scalar_real=wmax
    else if (name .eq. "wmin_local") then
      field_value%scalar_real=wmin
    else if (name .eq. "totqv_local") then
      field_value%scalar_real=totqv
    else if (name .eq. "totql_local") then
      field_value%scalar_real=totql
    else if (name .eq. "totqz_local") then
      field_value%scalar_real=totqz
    else if (name .eq. "qlmax_local") then
      field_value%scalar_real=qlmax
    else if (name .eq. "hqlmax_local") then
      field_value%scalar_int=hqlmax
    else if (name .eq. "cltop_local") then
      field_value%scalar_real=cltop
    else if (name .eq. "clbas_local") then
      field_value%scalar_real=clbas
    else if (name .eq. "cltop_av_local") then
      field_value%scalar_real=cltop_av
    else if (name .eq. "clbas_av_local") then
      field_value%scalar_real=clbas_av
    else if (name .eq. "ncl_col_local") then
      field_value%scalar_int=ncl_col
    end if
  end subroutine field_value_retrieval_callback
end module diagnostics_mod
