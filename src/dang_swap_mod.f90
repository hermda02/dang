module dang_swap_mod
    use healpix_types
    use pix_tools
    use fitstools
    use udgrade_nr
    use dang_util_mod
    use dang_param_mod
    use dang_linalg_mod
    use dang_component_mod
    use dang_data_mod
    implicit none


contains

  ! subroutine swap_bp_maps(dat,dpar,iteration,chain)
  !   type(params)                                    :: dpar
  !   type(data),                       intent(inout) :: dat
  !   integer(i4b),                     intent(in)    :: iteration
  !   character(len=512),               intent(in)    :: chain
  !   character(len=300), allocatable, dimension(:,:) :: bp_maps
  !   integer(i4b)                                    :: i, j
  !   character(len=6)                                :: iter_str
  !   real(dp), allocatable, dimension(:,:)           :: map, rms
    
    
  !   allocate(map(0:npix-1,3))
  !   allocate(rms(0:npix-1,3))
  !   allocate(bp_maps(dpar%numband,2))
    
  !   write(iter_str,'(i0.6)') iteration
    

  !   do j = 1, dpar%numband
  !      if (dpar%bp_map(j)) then
  !         bp_maps(j,1) = trim(dpar%bp_dir) // trim(dpar%band_label(j))//'_map_'//trim(chain)//'_n0064_60arcmin_k'//trim(iter_str)//'.fits'
  !         bp_maps(j,2) = trim(dpar%bp_dir) // trim(dpar%band_label(j))//'_rms_'//trim(chain)//'_n0064_60arcmin_k'//trim(iter_str)//'.fits'
  !         write(*,'(a,a,a)') 'Swapping band ', trim(dpar%band_label(j)), '.'
  !         !write(*,*) trim(bp_maps(j,1))                                                             
  !         !write(*,*) trim(bp_maps(j,2))                                                             
  !         call read_bintab(trim(bp_maps(j,1)),map,dat%npix,3,nullval,anynull,header=header)
  !         dat%sig_map(:,:,j) = map
  !         call read_bintab(trim(bp_maps(j,2)),rms,dat%npix,3,nullval,anynull,header=header)
  !         dat%rms_map(:,:,j) = rms
  !      end if
  !   end do


  subroutine swap_bp_maps(dat,dpar)
    type(params)                                    :: dpar
    type(data),                       intent(inout) :: dat
    character(len=512)                              :: chain_c
    character(len=300), allocatable, dimension(:,:) :: bp_maps
    integer(i4b)                                    :: i, j, iter_i, chain_i
    character(len=6)                                :: iter_str
    real(dp), allocatable, dimension(:,:)           :: map, rms
    real(dp)                                        :: norm
    double precision                                :: temp(1)

    allocate(map(0:npix-1,3))
    allocate(rms(0:npix-1,3))
    allocate(bp_maps(dpar%numband,2))

    do j = 1, dpar%numband
       if (dpar%bp_map(j)) then

          call RANDOM_SEED()
          call RANDOM_NUMBER(temp)
          
          norm    = temp(1)*dpar%num_chains
          chain_i = int(norm)+1
          chain_c = dpar%bp_chain_list(chain_i)
          
          call RANDOM_SEED()
          call RANDOM_NUMBER(temp)
          
          norm    = temp(1)*(dpar%bp_max-dpar%bp_burnin)
          iter_i  = int(norm)+1+dpar%bp_burnin
          
          write(iter_str,'(i0.6)') iter_i

          bp_maps(j,1) = trim(dpar%bp_dir) // trim(dpar%band_label(j))//'_map_'//trim(chain_c)//'_n0064_60arcmin_k'//trim(iter_str)  // '.fits'
          bp_maps(j,2) = trim(dpar%bp_dir) // trim(dpar%band_label(j))//'_rms_'//trim(chain_c)//'_n0064_60arcmin_k'//trim(iter_str) // '.fits'
          write(*,'(a,a,a)') 'Swapping band ', trim(dpar%band_label(j)), '.'
          call read_bintab(trim(bp_maps(j,1)),map,dat%npix,3,nullval,anynull,header=header)
          dat%sig_map(:,:,j) = map
          call read_bintab(trim(bp_maps(j,2)),rms,dat%npix,3,nullval,anynull,header=header)
          dat%rms_map(:,:,j) = rms
       end if
    end do

  end subroutine swap_bp_maps

end module dang_swap_mod
