module dang_swap_mod
    use healpix_types
    use pix_tools
    use fitstools
    use udgrade_nr
    use utility_mod
    use dang_param_mod
    use linalg_mod
    use dang_component_mod
    use dang_data_mod
    implicit none


contains


  subroutine swap_bp_maps(dat,param,iteration)
    type(params)                                    :: param
    type(data),                intent(inout)        :: dat
    integer(i4b),              intent(in)           :: iteration
    character(len=300), allocatable, dimension(:,:) :: bp_maps
    integer(i4b)                                    :: i, j
    character(len=6)                                :: iter_str
    real(dp), allocatable, dimension(:,:)           :: map, rms

    allocate(map(0:npix-1,3))
    allocate(rms(0:npix-1,3))
    allocate(bp_maps(param%numband,2))

    write(iter_str,'(i0.6)') iteration

    do j = 1, param%numband
       if (param%bp_map(j)) then
          bp_maps(j,1) = trim(param%bp_dir) // trim(param%dat_label(j))//'_map_n0064_60arcmin_k'//trim(iter_str)  // '.fits'
          bp_maps(j,2) = trim(param%bp_dir) // trim(param%dat_label(j))//'_rms_n0064_60arcmin_k'//trim(iter_str) // '.fits'
          write(*,'(a,a,a)') 'Swapping band ', trim(param%dat_label(j)), '.'
          call read_bintab(trim(bp_maps(j,1)),map,dat%npix,3,nullval,anynull,header=header)
          dat%sig_map(:,:,j) = map
          call read_bintab(trim(bp_maps(j,2)),rms,dat%npix,3,nullval,anynull,header=header)
          dat%rms_map(:,:,j) = rms
       end if
    end do

  end subroutine swap_bp_maps

end module dang_swap_mod
