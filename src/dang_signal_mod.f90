module dang_signal_mod
  use dang_component_mod
  use dang_param_mod
  use dang_cmb_comp_mod
  use dang_freefree_comp_mod
  use dang_hi_fit_comp_mod
  use dang_lognormal_comp_mod
  use dang_powlaw_comp_mod
  use dang_T_cmb_comp_mod
  use dang_mbb_comp_mod
  use dang_monopole_comp_mod
  use dang_template_comp_mod
  
contains
  
  subroutine initialize_components(dpar)
    implicit none
    type(dang_params)             :: dpar
    integer(i4b)                  :: i
    
    write(*,*) 'Initializing components'
    
    allocate(component_list(dpar%ncomp))
    
    do i = 1, dpar%ncomp
       if (trim(dpar%fg_type(i)) == 'mbb') then
          component_list(i)%p => dang_mbb_comp(dpar,i) 
       else if (trim(dpar%fg_type(i)) == 'power-law') then
          component_list(i)%p => dang_powlaw_comp(dpar,i)
       else if (trim(dpar%fg_type(i)) == 'freefree') then
          component_list(i)%p => dang_freefree_comp(dpar,i)
       else if (trim(dpar%fg_type(i)) == 'cmb') then
          component_list(i)%p => dang_cmb_comp(dpar,i)
       else if (trim(dpar%fg_type(i)) == 'T_cmb') then
          component_list(i)%p => dang_T_cmb_comp(dpar,i)
       else if (trim(dpar%fg_type(i)) == 'lognormal') then
          component_list(i)%p => dang_lognormal_comp(dpar,i)
       else if (trim(dpar%fg_type(i)) == 'template') then
          component_list(i)%p => dang_template_comp(dpar,i)
       else if (trim(dpar%fg_type(i)) == 'monopole') then
          component_list(i)%p => dang_monopole_comp(dpar,i)
       else if (trim(dpar%fg_type(i)) == 'hi_fit') then
          component_list(i)%p => dang_hi_fit_comp(dpar,i)
       end if
    end do

  end subroutine initialize_components

end module dang_signal_mod
