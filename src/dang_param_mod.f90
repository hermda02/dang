module dang_param_mod
  use healpix_types
  use dang_util_mod
  use hashtbl
  implicit none
  
  type, public :: dang_params
     ! Global parameters
     integer(i4b)                                  :: ngibbs           ! Number of Gibbs iterations
     integer(i4b)                                  :: nsample          ! For internal samplers (MH)
     integer(i4b)                                  :: iter_out         ! Out put maps every <- iterations
     integer(i4b)                                  :: cg_burnin        ! Number of burn in samples for cg chains
     integer(i4b)                                  :: cg_max           ! Maximum number of maps from cg chains
     integer(i4b)                                  :: num_chains       ! Number of cg chains used
     logical(lgt)                                  :: cg_swap          ! Do the cg map swapping?
     logical(lgt)                                  :: output_fg        ! Do we output the foregrounds at each frequency?
     character(len=512)                            :: outdir           ! Output directory
     character(len=512)                            :: cg_chains        ! cg chains, read as a string, moved to a 'list'
     character(len=16)                             :: ml_mode          ! 'sample' or 'optimize'
     character(len=16)                             :: solver           ! Linear system solver type 
     character(len=16)                             :: mode             ! 'dang' mode ('comp_sep', 'HI_fit')
     character(len=5)                              :: tqu              ! Which pol_type to sample
     real(dp)                                      :: cg_converge      ! cg convergence criterion 
     integer(i4b), allocatable, dimension(:)       :: pol_type         ! Points above to map number
     character(len=512), allocatable, dimension(:) :: cg_chain_list    ! This is where the cg_chains string goes
     
     ! Data parameters
     integer(i4b)                                    :: numband        ! Number of bands total in parameter file
     integer(i4b)                                    :: numinc         ! Number of bands to include in the fit
     character(len=512)                              :: datadir        ! Directory to look for bandfiles in
     character(len=512)                              :: cg_dir         ! Directory for cg swap maps
     character(len=512)                              :: mask_file      ! Mask filename
     character(len=512)                              :: offset_file    ! Offset init     
     character(len=512)                              :: gain_file      ! Gain init     
     character(len=512), allocatable, dimension(:)   :: band_calibrator! Band filename
     character(len=512), allocatable, dimension(:)   :: band_label     ! Band label
     character(len=512), allocatable, dimension(:)   :: band_mapfile   ! Band filename
     character(len=512), allocatable, dimension(:)   :: band_noisefile ! Band rms filename
     character(len=512), allocatable, dimension(:)   :: band_unit      ! Band units (uK_CMB, uK_RJ, MJy/sr)
     character(len=512), allocatable, dimension(:)   :: bp_id          ! Band bandpass type
     character(len=512), allocatable, dimension(:)   :: bp_file        ! Band bandpass filename
     real(dp),           allocatable, dimension(:)   :: band_nu        ! Band frequency (in GHz)
     real(dp),           allocatable, dimension(:)   :: init_gain      ! initial gain value for each band
     real(dp),           allocatable, dimension(:)   :: init_offset    ! initial offset value for each band
     logical(lgt),       allocatable, dimension(:)   :: cg_map         ! True false (know when to swap)
     logical(lgt),       allocatable, dimension(:)   :: fit_gain       ! Do we fit the gain for this band?
     logical(lgt),       allocatable, dimension(:)   :: fit_offs       ! Do we fit the offset for this band?
     logical(lgt),       allocatable, dimension(:)   :: band_inc       ! Is this band included?
     
     ! Component parameters
     integer(i4b)                                      :: ncomp          ! # of foregrounds
     integer(i4b)                                      :: ntemp          ! # of templates 
     
     character(len=512), allocatable, dimension(:)     :: temp_file      ! Template Filename
     character(len=100), allocatable, dimension(:)     :: temp_amps      ! Band filename
     character(len=512), allocatable, dimension(:)     :: temp_label     ! Template label
     character(len=10),  allocatable, dimension(:)     :: temp_polfit    ! Which poltypes are fit jointly for the template? ex. 'T', 'QU'
     logical(lgt),       allocatable, dimension(:,:)   :: temp_corr      ! Storing which bands should have templates fit
     ! integer(i4b),       allocatable, dimension(:)     :: temp_nfit      ! Number of bands fit for template i
     logical(lgt),       allocatable, dimension(:)     :: temp_sample    ! Storing which bands should have templates fit
     
     character(len=512), allocatable, dimension(:,:)   :: fg_amp_file    ! Fg amplitude input map
     logical(lgt),       allocatable, dimension(:)     :: fg_amp_samp    ! Logical - do we fit the amplitude for this guy?
     character(len=512), allocatable, dimension(:)     :: fg_filename    ! Init filename for diffuse components, template filename for templates
     integer(i4b),       allocatable, dimension(:)     :: fg_cg_group    ! Which cg group is this component in?
     real(dp),           allocatable, dimension(:,:)   :: fg_init        ! Initialized parameter value (fullsky)
     character(len=512), allocatable, dimension(:,:)   :: fg_ind_region  ! Fg spectral parameter input map
     character(len=512), allocatable, dimension(:,:)   :: fg_ind_lnl     ! Fg spectral parameter likelihood evaluation
     real(dp),           allocatable, dimension(:,:,:) :: fg_gauss       ! Fg gaussian sampling parameters
     character(len=512), allocatable, dimension(:)     :: fg_label       ! Fg label
     integer(i4b),       allocatable, dimension(:)     :: fg_nfit        ! How many bands are fit?
     real(dp),           allocatable, dimension(:)     :: fg_nu_ref      ! Fg reference frequency
     character(len=512), allocatable, dimension(:,:)   :: fg_prior_type  ! Fg spectral parameter input map
     logical(lgt),       allocatable, dimension(:,:)   :: fg_samp_spec   ! Logical - sample fg parameter?
     integer(i4b),       allocatable, dimension(:,:)   :: fg_samp_nside  ! Fg parameter nside sampling
     character(len=512), allocatable, dimension(:,:)   :: fg_spec_file   ! Fg spectral parameter input map
     character(len=10),  allocatable, dimension(:,:)   :: fg_spec_poltype ! Fg spectral parameter poltype
     real(dp),           allocatable, dimension(:,:)   :: fg_spec_step   ! Fg index step-size of MH
     logical(lgt),       allocatable, dimension(:,:)   :: fg_spec_tune   ! Logical - tune spectral parameter step size?
     logical(lgt),       allocatable, dimension(:,:)   :: fg_temp_corr   ! Logical - do we template correct this band?
     character(len=512), allocatable, dimension(:)     :: fg_type        ! Fg type (power-law feks)
     real(dp),           allocatable, dimension(:,:,:) :: fg_uni         ! Fg sampling bounds
     
     logical(lgt)                                      :: joint_sample   ! Logical - jointly sample fg amplitudes
     logical(lgt)                                      :: joint_pol      ! Points to which Stokes are jointly sampled
     character(len=512), allocatable, dimension(:)     :: joint_comp     ! Joint sampler components
     
     real(dp)                                          :: thresh         ! Threshold for the HI fitting (sample pixels under thresh)
     character(len=512)                                :: HI_file        ! HI map filename
     character(len=512)                                :: HI_Td_init     ! HI fitting dust temp estimate
     real(dp)                                          :: HI_Td_mean     ! HI Temperature prior mean
     real(dp)                                          :: HI_Td_std      ! HI Temperature prior std
     real(dp)                                          :: HI_Td_step     ! Td sampling step size

     ! cg group parameters
     integer(i4b)                                      :: ncggroup       ! How many CG groups do we expect to read?
     logical(lgt),       allocatable, dimension(:)     :: cg_group_sample ! Do we sample that cg group?
     integer(i4b),       allocatable, dimension(:)     :: cg_max_iter    ! What is the maximum # of CG iterations we perform?
     real(dp),           allocatable, dimension(:)     :: cg_convergence ! What is the convergence criterion for the CG group?
     character(len=10),  allocatable, dimension(:)     :: cg_poltype     ! What is the CG group poltype?

  end type dang_params
  
contains

  subroutine get_file_length(filename,length)
    implicit none
    character(len=512), intent(in)  :: filename
    integer(i4b),       intent(out) :: length
    integer(i4b), parameter    :: maxdepth = 256
    integer(i4b)               :: depth, units(maxdepth), i
    character(len=512)         :: key, value, filenames(maxdepth), line
    
    length = 0
    depth = 1
    units(depth) = getlun()
    filenames(depth) = filename
    open(units(depth),file=trim(filename),status="old")!,err=4)
    do while(depth >= 1)
       read(units(depth),*,end=1) key
       if (key(1:1)=='#') cycle
       backspace(units(depth))
       
       if (key(1:1)=='@') then
          if(key == '@INCLUDE') then
             ! Recurse into the new file                                                             
             read(units(depth),*,end=1) key, value
             depth=depth+1
             units(depth) = getlun()
             filenames(depth) = value
             open(units(depth),file=value,status="old")!,err=2)
          else
             stop
          end if
       else
          read(units(depth),fmt="(a)") line
          !if we get here we have read a new line from the parameter file(s)              
          length = length + 1
       end if
       cycle
       ! We get here if we reached the end of a file. Close it and                                   
       ! return to the file above.                                                                   
1      close(units(depth))
       depth = depth-1
    end do
    return
    
  end subroutine get_file_length
  
  subroutine read_param_file(par)
    implicit none
    type(hash_tbl_sll)                            :: htable
    type(dang_params), intent(inout)              :: par
    integer(i4b)                                  :: parfile_len, i
    character(len=512)                            :: paramfile
    character(len=512), allocatable, dimension(:) :: parfile_cache
    
    call getarg(1,paramfile)

    write(*,*) ''
    if (rank == master) write(*,*) 'Reading parameters from ', trim(paramfile)
    write(*,*) ''
    write(*,*) '--------------------------'
   
    call get_file_length(paramfile,parfile_len)
    allocate(parfile_cache(parfile_len))
    call read_paramfile_to_ascii(paramfile,parfile_cache)
    
    !Initialize a hash table                                                                         
    call init_hash_tbl_sll(htable,tbl_len=10*parfile_len)
    ! Put the parameter file into the hash table                                                     
    call put_ascii_into_hashtable(parfile_cache,htable)
    
    call read_global_params(htable,par)    
    call read_data_params(htable,par)
    call read_comp_params(htable,par)
    deallocate(parfile_cache)
    write(*,*) '--------------------------'

  end subroutine read_param_file
  
  subroutine read_paramfile_to_ascii(paramfile,paramfile_cache)
    implicit none
    character(len=512),                            intent(in)    :: paramfile
    character(len=512), allocatable, dimension(:), intent(inout) :: paramfile_cache
    integer(i4b), parameter   :: maxdepth = 256
    integer(i4b)              :: depth, units(maxdepth), line_nr, paramfile_len, i
    character(len=512)        :: key, value, filenames(maxdepth), line
    
    line_nr = 0
    depth   = 1
    units(depth) = getlun()
    filenames(depth) = paramfile
    open(units(depth),file=trim(paramfile),status="old")!,err=4)
    do while(depth >= 1)
       read(units(depth),*,end=1) key
       if (key(1:1)=='#') cycle
       backspace(units(depth))
       
       if (key(1:1)=='@') then
          if(key == '@INCLUDE') then
             ! Recurse into the new file                                                             
             read(units(depth),*,end=1) key, value
             depth=depth+1
             units(depth) = getlun()
             filenames(depth) = value
             open(units(depth),file=value,status="old")!,err=2)
          else
             stop
          end if
       else
          read(units(depth),fmt="(a)") line
          !if we get here we have read a new line from the parameter file(s)                         
          line_nr = line_nr + 1
          write(paramfile_cache(line_nr),fmt="(a)") line
       end if
       cycle
1      close(units(depth))
       depth = depth-1
    end do
    return
    
  end subroutine read_paramfile_to_ascii
  
  subroutine put_ascii_into_hashtable(asciitbl,htbl)
    implicit none
    character(len=512), allocatable, dimension(:), intent(in) :: asciitbl
    type(hash_tbl_sll), intent(inout) :: htbl
    character(len=512) :: key, val
    character(len=256) :: toks(2)
    integer            :: i, n
    do i = 1,size(asciitbl)
       call get_tokens(trim(asciitbl(i)), "=", group="''" // '""', maxnum=2, toks=toks, num=n)
       if(n < 2) then ! only need the lines where one has 'key'='value'                                                                   
          cycle
       end if
       key = get_token(toks(1), " ", 1, group="''" // '""')
       val = get_token(toks(2), " ", 1, group="''" // '""')
       call tolower(key)  ! we don't differentiate btw. upper and lower case                                                              
       if (key=="") cycle ! we don't need blank lines                                                                                      
       call put_hash_tbl_sll(htbl,trim(key),trim(val))
    end do
    return
    
    write(*,*) "Error: Cannot read ascii line:", i, "line = '" // trim(asciitbl(i)) // "'"
    stop
    
  end subroutine put_ascii_into_hashtable
  
  ! read parameter from input argument or hash table                                                 
  subroutine get_parameter_hashtable(htbl, parname, len_itext, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    type(hash_tbl_sll), intent(in) :: htbl 
    character(len=*),   intent(in) :: parname
    integer(i4b),         optional :: len_itext
    integer(i4b),         optional :: par_int
    character(len=*),     optional :: par_char
    character(len=*),     optional :: par_string
    real(sp),             optional :: par_sp
    real(dp),             optional :: par_dp
    logical(lgt),         optional :: par_lgt
    logical(lgt),         optional :: par_present
    character(len=*),     optional :: desc
    
    logical(lgt)               :: found
    
    call get_parameter_from_hash(htbl, parname, len_itext, par_int, &
         & par_char, par_string, par_sp, par_dp, par_lgt, par_present, desc)
  end subroutine get_parameter_hashtable
  
  ! getting parameter value from hash table                                                          
  subroutine get_parameter_from_hash(htbl, parname, len_itext, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    type(hash_tbl_sll), intent(in) :: htbl
    character(len=*),   intent(in) :: parname
    integer(i4b),     optional :: len_itext
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc
    character(len=256)         :: key
    character(len=:), ALLOCATABLE   :: itext,jtext
    CHARACTER(len=:), ALLOCATABLE   :: val,val2,val3
    integer(i4b)                    :: i,j
    
    key=adjustl(trim(parname))
    call tolower(key)
    call get_hash_tbl_sll(htbl,trim(key),val)
    if (.not. allocated(val)) then
       goto 1
       if (.not. present(len_itext)) goto 1
       allocate(character(len=len_itext) :: itext,jtext)
       itext=key(len(trim(key))-(len_itext-1):len(trim(key)))
       call get_hash_tbl_sll(htbl,'band_default_params'//trim(itext),val2)
       if (allocated(val2)) then
          read(val2,*) j
          if (j /= 0) then
             call int2string(j, jtext)
             call get_hash_tbl_sll(htbl,'band_default_params'//trim(jtext),val3)
             if (allocated(val3)) then
                read(val3,*) i
                if (i /= 0) goto 2
             end if
             call get_hash_tbl_sll(htbl,key(1:len(trim(key))-len_itext)//trim(jtext),val)
             if (.not. allocated(val)) goto 3
          else
             goto 1
          end if
       else
          goto 1
       end if
       deallocate(itext,jtext)
    end if
    
    if (present(par_int)) then
       read(val,*) par_int
    elseif (present(par_char)) then
       read(val,*) par_char
    elseif (present(par_string)) then
       read(val,*) par_string
    elseif (present(par_sp)) then
       read(val,*) par_sp
    elseif (present(par_dp)) then
       read(val,*) par_dp
    elseif (present(par_lgt)) then
       read(val,*) par_lgt
    else
       write(*,*) "get_parameter: Reached unreachable point!"
    end if
    
    deallocate(val)
    return
    
1   write(*,*) "Error: Could not find parameter '" // trim(parname) // "'"
    write(*,*) ""
    stop
    
    
2   write(*,*) "Error: Recursive default parameters, bands " // &
         & trim(jtext) // " and " //trim(itext)
    write(*,*) ""
    stop
    
3   write(*,*) "Error: Could not find parameter '" // trim(parname) // &
         & "' from default '"//key(1:len(trim(key))-len_itext)//trim(jtext)//"'"
    write(*,*) ""
    stop
    
  end subroutine get_parameter_from_hash
  
  subroutine read_global_params(htbl,par)
    implicit none
    
    type(hash_tbl_sll), intent(in)    :: htbl
    type(dang_params),       intent(inout) :: par
    integer(i4b)                      :: pol_count
    
    integer(i4b)     :: i, j, n, len_itext
    character(len=2) :: itext
    character(len=2) :: jtext
    
    write(*,*) "Read global parameters."
    
    call get_parameter_hashtable(htbl, 'OUTPUT_DIRECTORY', par_string=par%outdir)
    call get_parameter_hashtable(htbl, 'NUMGIBBS', par_int=par%ngibbs)
    call get_parameter_hashtable(htbl, 'NUMSAMPLE', par_int=par%nsample)
    call get_parameter_hashtable(htbl, 'OUTPUT_ITER', par_int=par%iter_out)
    call get_parameter_hashtable(htbl, 'OUTPUT_COMPS', par_lgt=par%output_fg)
    call get_parameter_hashtable(htbl, 'SOLVER_TYPE', par_string=par%solver)
    call get_parameter_hashtable(htbl, 'ML_MODE', par_string=par%ml_mode)
    call get_parameter_hashtable(htbl, 'TQU', par_string=par%tqu)
    call get_parameter_hashtable(htbl, 'CG_SWAP',par_lgt=par%cg_swap)

    ! Gather parameters for conjugate gradient solver
    call get_parameter_hashtable(htbl, 'NUM_CG_GROUPS', par_int=par%ncggroup) 
    n = par%ncggroup

    allocate(par%cg_group_sample(n))
    allocate(par%cg_max_iter(n))
    allocate(par%cg_convergence(n))
    allocate(par%cg_poltype(n))

    ! Load the CG group specific parameters
    do i = 1, n
       call int2string(i, itext)
       call get_parameter_hashtable(htbl, 'CG_GROUP_SAMPLE'//itext, len_itext=len_itext, &
            par_lgt=par%cg_group_sample(i))
       call get_parameter_hashtable(htbl, 'CG_GROUP_MAX_ITER'//itext, len_itext=len_itext, &
            par_int=par%cg_max_iter(i))
       call get_parameter_hashtable(htbl, 'CG_CONVERGE_THRESH'//itext, len_itext=len_itext, &
            par_dp=par%cg_convergence(i))
       call get_parameter_hashtable(htbl, 'CG_POLTYPE'//itext, len_itext=len_itext, &
            par_string=par%cg_poltype(i))
    end do

    write(*,*) size(par%cg_poltype), par%cg_convergence(1), par%cg_poltype(1)
    

    if (par%cg_swap) then
         call get_parameter_hashtable(htbl, 'CG_BURN_IN',par_int=par%cg_burnin)
         call get_parameter_hashtable(htbl, 'CG_MAX_ITER',par_int=par%cg_max)
         call get_parameter_hashtable(htbl, 'CG_DIRECTORY',par_string=par%cg_dir)
         call get_parameter_hashtable(htbl, 'CG_CHAINS_LIST',par_string=par%cg_chains)
         call get_parameter_hashtable(htbl, 'CG_NUM_CHAINS',par_int=par%num_chains)
    end if

    if (par%num_chains /= 0) then
       allocate(par%cg_chain_list(par%num_chains))
    
       call delimit_string(par%cg_chains,',',par%cg_chain_list)
    end if
    
    ! Surely an inefficient way to decide which maps to use (T -> 1, Q -> 2, U -> 3), but it works
    pol_count = 0
    if (index(par%tqu,'T') /= 0) then
       pol_count = pol_count + 1
    end if
    if (index(par%tqu,'Q') /= 0) then
       pol_count = pol_count + 1
    end if
    if (index(par%tqu,'U') /= 0) then
       pol_count = pol_count + 1
    end if
    allocate(par%pol_type(pol_count))
    
    pol_count = 0
    if (index(par%tqu,'T') /= 0) then
       pol_count = pol_count + 1
       par%pol_type(pol_count) = 1
    end if
    if (index(par%tqu,'Q') /= 0) then
       pol_count = pol_count + 1
       par%pol_type(pol_count) = 2
    end if
    if (index(par%tqu,'U') /= 0) then
       pol_count = pol_count + 1
       par%pol_type(pol_count) = 3
    end if
    
    par%outdir = trim(par%outdir) // '/'
    
    write(*,*) size(par%cg_poltype), par%cg_convergence(1), par%cg_poltype(1)
    
  end subroutine read_global_params
  
  subroutine read_data_params(htbl,par)
    implicit none
    
    type(hash_tbl_sll), intent(in)    :: htbl
    type(dang_params),       intent(inout) :: par
    
    integer(i4b)     :: i, j, n, n2, len_itext
    character(len=3) :: itext
    character(len=2) :: jtext
    
    write(*,*) "Read data parameters."
    
    len_itext = len(trim(itext))
    
    call get_parameter_hashtable(htbl, 'NUMBAND',    par_int=par%numband)
    call get_parameter_hashtable(htbl, 'NUMINCLUDE',    par_int=par%numinc)
    call get_parameter_hashtable(htbl, 'DATA_DIRECTORY', par_string=par%datadir)
    call get_parameter_hashtable(htbl, 'MASKFILE', par_string=par%mask_file)

    call get_parameter_hashtable(htbl, 'BAND_OFFSET_FILE', par_string=par%offset_file)
    call get_parameter_hashtable(htbl, 'BAND_GAIN_FILE', par_string=par%gain_file)
    
    n  = par%numband
    n2 = par%numinc
    
    allocate(par%band_inc(n))
    allocate(par%band_mapfile(n2),par%band_label(n2))
    allocate(par%band_noisefile(n2),par%band_nu(n2))
    allocate(par%band_calibrator(n2))
    allocate(par%bp_id(n2),par%bp_file(n2))
    allocate(par%band_unit(n2))
    allocate(par%cg_map(n2))
    
    allocate(par%init_gain(n2))
    allocate(par%init_offset(n2))
    allocate(par%fit_gain(n2))
    allocate(par%fit_offs(n2))
    
    ! Set up this way so as to only load information about included bands
    j = 0
    do i = 1, n
       call int2string(i, itext)
       call get_parameter_hashtable(htbl, 'INCLUDE_BAND'//itext, len_itext=len_itext, par_lgt=par%band_inc(i))
       ! Cycle past any bands not included
       if (.not. par%band_inc(i)) then
          cycle
       else
          j = j + 1
       end if
       call get_parameter_hashtable(htbl, 'BP_TYPE' // itext, par_string=par%bp_id(j))
       if (trim(par%bp_id(j)) == 'delta') then
          par%bp_file(j) = ''
       else
          call get_parameter_hashtable(htbl, 'BP_FILE' // itext, par_string=par%bp_file(j))
       end if

       ! call get_parameter_hashtable(htbl, 'BAND_CALIBRATOR'//itext, len_itext=len_itext, par_string=par%band_calibrator(j))
       call get_parameter_hashtable(htbl, 'BAND_LABEL'//itext, len_itext=len_itext, par_string=par%band_label(j))
       call get_parameter_hashtable(htbl, 'BAND_FILE'//itext, len_itext=len_itext, par_string=par%band_mapfile(j))
       call get_parameter_hashtable(htbl, 'BAND_RMS'//itext, len_itext=len_itext, par_string=par%band_noisefile(j))
       call get_parameter_hashtable(htbl, 'BAND_FREQ'//itext, len_itext=len_itext, par_dp=par%band_nu(j))
       call get_parameter_hashtable(htbl, 'BAND_UNIT'//itext, len_itext=len_itext, par_string=par%band_unit(j))
       call get_parameter_hashtable(htbl, 'BAND_FIT_GAIN'//itext, len_itext=len_itext, par_lgt=par%fit_gain(j))
       ! call get_parameter_hashtable(htbl, 'BAND_FIT_OFFSET'//itext, len_itext=len_itext, par_lgt=par%fit_offs(j))
       call get_parameter_hashtable(htbl, 'BAND_CG'//itext, len_itext=len_itext, par_lgt=par%cg_map(j))
    end do
  end subroutine read_data_params
  
  subroutine read_comp_params(htbl,par)
    implicit none
    
    type(hash_tbl_sll), intent(in)    :: htbl
    type(dang_params),       intent(inout) :: par
    
    integer(i4b)     :: i, j, n, n2, n3, n4
    integer(i4b)     :: len_itext, len_jtext
    character(len=2) :: itext
    character(len=3) :: jtext
    
    write(*,*) "Read component parameters."
    
    len_itext = len(trim(itext))
    len_jtext = len(trim(jtext))

    call get_parameter_hashtable(htbl, 'NUMCOMPS', par_int=par%ncomp)
    call get_parameter_hashtable(htbl, 'NUMTEMPS', par_int=par%ntemp)
   
    n  = par%ncomp
    
    allocate(par%fg_amp_file(n,1))
    allocate(par%fg_amp_samp(n))
    allocate(par%fg_cg_group(n))
    allocate(par%fg_filename(n))
    allocate(par%fg_gauss(n,2,2),par%fg_uni(n,2,2))
    allocate(par%fg_ind_region(n,2))
    allocate(par%fg_ind_lnl(n,2))
    allocate(par%fg_init(n,2))
    allocate(par%fg_label(n),par%fg_type(n),par%fg_nu_ref(n))
    allocate(par%fg_nfit(n))
    allocate(par%fg_prior_type(n,2))
    allocate(par%fg_samp_nside(n,2),par%fg_samp_spec(n,2))
    allocate(par%fg_spec_file(n,2))
    allocate(par%fg_spec_poltype(n,2))
    allocate(par%fg_spec_step(n,2))
    allocate(par%fg_spec_tune(n,2))
    allocate(par%fg_temp_corr(n,par%numband))
        
    allocate(par%temp_file(par%ntemp))
    allocate(par%temp_amps(par%ntemp))
    
    ! par%temp_nfit = 0
    par%fg_nfit   = 0
    
    ! Loop over and load all component information
    do i = 1, n
       call int2string(i, itext)
       call get_parameter_hashtable(htbl, 'COMP_LABEL'//itext, len_itext=len_itext, par_string=par%fg_label(i))
       call get_parameter_hashtable(htbl, 'COMP_TYPE'//itext, len_itext=len_itext, par_string=par%fg_type(i)) 
       call get_parameter_hashtable(htbl, 'COMP_CG_GROUP'//itext, len_itext=len_itext, par_int=par%fg_cg_group(i))
       if (trim(par%fg_type(i)) /= 'monopole') then
          call get_parameter_hashtable(htbl, 'COMP_FILENAME'//itext, len_itext=len_itext, par_string=par%fg_filename(i))
       end if
       if (trim(par%fg_type(i)) /= 'template' .and. trim(par%fg_type(i)) /= 'hi_fit' .and. trim(par%fg_type(i)) /= 'monopole') then
          call get_parameter_hashtable(htbl, 'COMP_AMP_SAMPLE'//itext, len_itext=len_itext, par_lgt=par%fg_amp_samp(i))
          call get_parameter_hashtable(htbl, 'COMP_REF_FREQ'//itext, len_itext=len_itext, par_dp=par%fg_nu_ref(i))
          if (par%fg_nu_ref(i) < 1d7) then
             par%fg_nu_ref(i) = par%fg_nu_ref(i)*1d9
          end if
       end if
       if (trim(par%fg_type(i)) == 'power-law') then
          call read_power_law(par,htbl,i)
       else if (trim(par%fg_type(i)) == 'cmb') then
          cycle
       else if (trim(par%fg_type(i)) == 'T_cmb') then
          call read_T_cmb(par,htbl,i)
       else if (trim(par%fg_type(i)) == 'freefree') then
          call read_freefree(par,htbl,i)
       else if (trim(par%fg_type(i)) == 'lognormal') then
          call read_lognormal(par,htbl,i)
       else if (trim(par%fg_type(i)) == 'mbb') then
          call read_mbb(par,htbl,i)
       else if (trim(par%fg_type(i)) == 'template') then
          call read_template(par,htbl,i)
       else if (trim(par%fg_type(i)) == 'monopole') then
          call read_monopole(par,htbl,i)
       else if (trim(par%fg_type(i)) == 'hi_fit') then
          call read_hi_fit(par,htbl,i)
       end if
    end do
    
  end subroutine read_comp_params

  subroutine read_T_cmb(par,htbl,comp)
    implicit none
    type(hash_tbl_sll),  intent(in)    :: htbl
    type(dang_params),   intent(inout) :: par
    integer(i4b),        intent(in)    :: comp

    integer(i4b)     :: len_itext
    character(len=2) :: itext
    
    len_itext = len(trim(itext))
    
    call int2string(comp, itext)

    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR_GAUSS_MEAN'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,1,1))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR_GAUSS_STD'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,1,2))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR_UNI_LOW'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,1,1))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR_UNI_HIGH'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,1,2))
    call get_parameter_hashtable(htbl, 'COMP_T_POLTYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_poltype(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T'//itext, len_itext=len_itext, par_dp=par%fg_init(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_SAMP_NSIDE'//itext, len_itext=len_itext,&
         par_int=par%fg_samp_nside(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_SAMPLE'//itext, len_itext=len_itext,&
         par_lgt=par%fg_samp_spec(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_INPUT_MAP'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_file(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_REGION'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_region(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_LNL_TYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_lnl(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR'//itext, len_itext=len_itext,&
         par_string=par%fg_prior_type(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_STEPSIZE'//itext,len_itext=len_itext,&
         par_dp=par%fg_spec_step(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_TUNE_STEPSIZE'//itext,len_itext=len_itext,&
         par_lgt=par%fg_spec_tune(comp,1))

  end subroutine read_T_cmb

  subroutine read_power_law(par,htbl,comp)
    implicit none
    type(hash_tbl_sll),  intent(in)    :: htbl
    type(dang_params),   intent(inout) :: par
    integer(i4b),        intent(in)    :: comp

    integer(i4b)     :: len_itext
    character(len=2) :: itext
    
    len_itext = len(trim(itext))
    
    call int2string(comp, itext)

    call get_parameter_hashtable(htbl, 'COMP_BETA_PRIOR_GAUSS_MEAN'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,1,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_PRIOR_GAUSS_STD'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,1,2))
    call get_parameter_hashtable(htbl, 'COMP_BETA_PRIOR_UNI_LOW'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,1,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_PRIOR_UNI_HIGH'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,1,2))
    call get_parameter_hashtable(htbl, 'COMP_BETA'//itext, len_itext=len_itext, par_dp=par%fg_init(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_SAMP_NSIDE'//itext, len_itext=len_itext,&
         par_int=par%fg_samp_nside(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_SAMPLE'//itext, len_itext=len_itext,&
         par_lgt=par%fg_samp_spec(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_POLTYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_poltype(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_INPUT_MAP'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_file(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_REGION'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_region(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_LNL_TYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_lnl(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_PRIOR'//itext, len_itext=len_itext,&
         par_string=par%fg_prior_type(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_STEPSIZE'//itext,len_itext=len_itext,&
         par_dp=par%fg_spec_step(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_TUNE_STEPSIZE'//itext,len_itext=len_itext,&
         par_lgt=par%fg_spec_tune(comp,1))

  end subroutine read_power_law

  subroutine read_freefree(par,htbl,comp)
    implicit none
    type(hash_tbl_sll),  intent(in)    :: htbl
    type(dang_params),   intent(inout) :: par
    integer(i4b),        intent(in)    :: comp

    integer(i4b)     :: len_itext
    character(len=2) :: itext
    
    len_itext = len(trim(itext))
    
    call int2string(comp, itext)

    call get_parameter_hashtable(htbl, 'COMP_T_E_PRIOR_GAUSS_MEAN'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,1,1))
    call get_parameter_hashtable(htbl, 'COMP_T_E_PRIOR_GAUSS_STD'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,1,2))
    call get_parameter_hashtable(htbl, 'COMP_T_E_PRIOR_UNI_LOW'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,1,1))
    call get_parameter_hashtable(htbl, 'COMP_T_E_PRIOR_UNI_HIGH'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,1,2))
    call get_parameter_hashtable(htbl, 'COMP_T_E'//itext, len_itext=len_itext, par_dp=par%fg_init(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_E_SAMP_NSIDE'//itext, len_itext=len_itext,&
         par_int=par%fg_samp_nside(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_E_SAMPLE'//itext, len_itext=len_itext,&
         par_lgt=par%fg_samp_spec(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_E_POLTYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_poltype(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_E_INPUT_MAP'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_file(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_E_REGION'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_region(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_E_LNL_TYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_lnl(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_E_PRIOR'//itext, len_itext=len_itext,&
         par_string=par%fg_prior_type(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_E_STEPSIZE'//itext,len_itext=len_itext,&
         par_dp=par%fg_spec_step(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_E_TUNE_STEPSIZE'//itext,len_itext=len_itext,&
         par_lgt=par%fg_spec_tune(comp,1))

  end subroutine read_freefree

  subroutine read_template(par,htbl,comp)
    implicit none
    type(hash_tbl_sll),  intent(in)    :: htbl
    type(dang_params),   intent(inout) :: par
    integer(i4b),        intent(in)    :: comp

    integer(i4b)     :: len_itext, len_jtext, i, j
    character(len=2) :: itext
    character(len=3) :: jtext
    
    len_itext = len(trim(itext))
    len_jtext = len(trim(jtext))
    
    call int2string(comp, itext)

    ! do i = 1, par%ntemp
    call get_parameter_hashtable(htbl, 'COMP_FILENAME'//itext, par_string=par%temp_file(comp))
    call get_parameter_hashtable(htbl, 'COMP_AMP_FILE'//itext, par_string=par%temp_amps(comp))
    call get_parameter_hashtable(htbl, 'COMP_POLTYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_poltype(comp,1))
    ! end do
    do j = 1, par%numband
       call int2string(j,jtext)
       call get_parameter_hashtable(htbl, 'COMP'//trim(itext)//'_FIT'//jtext,&
            len_itext=len_jtext,par_lgt=par%fg_temp_corr(comp,j))
       if (par%fg_temp_corr(comp,j)) then
          par%fg_nfit(comp) = par%fg_nfit(comp) + 1
       end if
    end do

  end subroutine read_template

  subroutine read_monopole(par,htbl,comp)
    implicit none
    type(hash_tbl_sll),  intent(in)    :: htbl
    type(dang_params),   intent(inout) :: par
    integer(i4b),        intent(in)    :: comp

    integer(i4b)     :: len_itext, len_jtext, i, j
    character(len=2) :: itext
    character(len=3) :: jtext
    
    len_itext = len(trim(itext))
    len_jtext = len(trim(jtext))
    
    call int2string(comp, itext)

    do j = 1, par%numband
       call int2string(j,jtext)
       call get_parameter_hashtable(htbl, 'COMP'//trim(itext)//'_FIT'//jtext,&
            len_itext=len_jtext,par_lgt=par%fg_temp_corr(comp,j))
       if (par%fg_temp_corr(comp,j)) then
          par%fg_nfit(comp) = par%fg_nfit(comp) + 1
       end if
    end do

  end subroutine read_monopole


  subroutine read_mbb(par,htbl,comp)
    implicit none
    type(hash_tbl_sll),  intent(in)    :: htbl
    type(dang_params),   intent(inout) :: par
    integer(i4b),        intent(in)    :: comp

    integer(i4b)     :: len_itext
    character(len=2) :: itext
    
    len_itext = len(trim(itext))
    
    call int2string(comp, itext)
    
    call get_parameter_hashtable(htbl, 'COMP_BETA_PRIOR_GAUSS_MEAN'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,1,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_PRIOR_GAUSS_STD'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,1,2))
    call get_parameter_hashtable(htbl, 'COMP_BETA_PRIOR_UNI_LOW'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,1,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_PRIOR_UNI_HIGH'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,1,2))
    call get_parameter_hashtable(htbl, 'COMP_BETA_POLTYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_poltype(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA'//itext, len_itext=len_itext, par_dp=par%fg_init(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_INPUT_MAP'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_file(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_SAMP_NSIDE'//itext, len_itext=len_itext,&
         par_int=par%fg_samp_nside(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_SAMPLE'//itext, len_itext=len_itext,&
         par_lgt=par%fg_samp_spec(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_REGION'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_region(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_LNL_TYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_lnl(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_PRIOR'//itext, len_itext=len_itext,&
         par_string=par%fg_prior_type(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_STEPSIZE'//itext,len_itext=len_itext,&
         par_dp=par%fg_spec_step(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_BETA_TUNE_STEPSIZE'//itext,len_itext=len_itext,&
         par_lgt=par%fg_spec_tune(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR_GAUSS_MEAN'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,2,1))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR_GAUSS_STD'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,2,2))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR_UNI_LOW'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,2,1))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR_UNI_HIGH'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,2,2))
    call get_parameter_hashtable(htbl, 'COMP_T_POLTYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_poltype(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_T'//itext, len_itext=len_itext, par_dp=par%fg_init(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_T_SAMP_NSIDE'//itext, len_itext=len_itext,&
         par_int=par%fg_samp_nside(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_T_SAMPLE'//itext, len_itext=len_itext,&
         par_lgt=par%fg_samp_spec(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_T_INPUT_MAP'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_file(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_T_REGION'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_region(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_T_LNL_TYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_lnl(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR'//itext, len_itext=len_itext,&
         par_string=par%fg_prior_type(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_T_STEPSIZE'//itext,len_itext=len_itext,&
         par_dp=par%fg_spec_step(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_T_TUNE_STEPSIZE'//itext,len_itext=len_itext,&
         par_lgt=par%fg_spec_tune(comp,2))

  end subroutine read_mbb

  subroutine read_lognormal(par,htbl,comp)
    implicit none
    type(hash_tbl_sll),  intent(in)    :: htbl
    type(dang_params),   intent(inout) :: par
    integer(i4b),        intent(in)    :: comp

    integer(i4b)     :: len_itext
    character(len=2) :: itext
    
    len_itext = len(trim(itext))
    
    call int2string(comp, itext)
    
    call get_parameter_hashtable(htbl, 'COMP_NU_P_PRIOR_GAUSS_MEAN'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,1,1))
    call get_parameter_hashtable(htbl, 'COMP_NU_P_PRIOR_GAUSS_STD'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,1,2))
    call get_parameter_hashtable(htbl, 'COMP_NU_P_PRIOR_UNI_LOW'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,1,1))
    call get_parameter_hashtable(htbl, 'COMP_NU_P_PRIOR_UNI_HIGH'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,1,2))
    call get_parameter_hashtable(htbl, 'COMP_NU_P_POLTYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_poltype(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_NU_P'//itext, len_itext=len_itext, par_dp=par%fg_init(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_NU_P_INPUT_MAP'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_file(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_NU_P_SAMP_NSIDE'//itext, len_itext=len_itext,&
         par_int=par%fg_samp_nside(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_NU_P_SAMPLE'//itext, len_itext=len_itext,&
         par_lgt=par%fg_samp_spec(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_NU_P_REGION'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_region(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_NU_P_LNL_TYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_lnl(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_NU_P_PRIOR'//itext, len_itext=len_itext,&
         par_string=par%fg_prior_type(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_NU_P_STEPSIZE'//itext,len_itext=len_itext,&
         par_dp=par%fg_spec_step(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_NU_P_TUNE_STEPSIZE'//itext,len_itext=len_itext,&
         par_lgt=par%fg_spec_tune(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_PRIOR_GAUSS_MEAN'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,2,1))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_PRIOR_GAUSS_STD'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,2,2))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_PRIOR_UNI_LOW'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,2,1))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_PRIOR_UNI_HIGH'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,2,2))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_POLTYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_poltype(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_W_AME'//itext, len_itext=len_itext, par_dp=par%fg_init(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_SAMP_NSIDE'//itext, len_itext=len_itext,&
         par_int=par%fg_samp_nside(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_SAMPLE'//itext, len_itext=len_itext,&
         par_lgt=par%fg_samp_spec(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_INPUT_MAP'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_file(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_REGION'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_region(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_LNL_TYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_lnl(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_PRIOR'//itext, len_itext=len_itext,&
         par_string=par%fg_prior_type(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_STEPSIZE'//itext,len_itext=len_itext,&
         par_dp=par%fg_spec_step(comp,2))
    call get_parameter_hashtable(htbl, 'COMP_W_AME_TUNE_STEPSIZE'//itext,len_itext=len_itext,&
         par_lgt=par%fg_spec_tune(comp,2))

  end subroutine read_lognormal

  subroutine read_hi_fit(par,htbl,comp)
    implicit none
    type(hash_tbl_sll),  intent(in)    :: htbl
    type(dang_params),   intent(inout) :: par
    integer(i4b),        intent(in)    :: comp

    integer(i4b)     :: len_itext, len_jtext, j
    character(len=2) :: itext
    character(len=3) :: jtext
    
    len_itext = len(trim(itext))
    len_jtext = len(trim(jtext))

    call get_parameter_hashtable(htbl, 'HI_THRESH', par_dp=par%thresh)
    call get_parameter_hashtable(htbl, 'HI_FILE',   par_string=par%HI_file)
    
    call int2string(comp, itext)

    call get_parameter_hashtable(htbl, 'COMP_LABEL'//itext, len_itext=len_itext, par_string=par%fg_label(comp))
    call get_parameter_hashtable(htbl, 'COMP_TYPE'//itext, len_itext=len_itext, par_string=par%fg_type(comp))
    call get_parameter_hashtable(htbl, 'COMP_REF_FREQ'//itext, len_itext=len_itext, par_dp=par%fg_nu_ref(comp))
    if (par%fg_nu_ref(comp) < 1d7) then
       par%fg_nu_ref(comp) = par%fg_nu_ref(comp)*1d9
    end if
    call get_parameter_hashtable(htbl, 'COMP_CG_GROUP'//itext, len_itext=len_itext, par_int=par%fg_cg_group(comp))
    call get_parameter_hashtable(htbl, 'COMP_FILENAME'//itext, len_itext=len_itext, par_string=par%fg_filename(comp))
    call get_parameter_hashtable(htbl, 'COMP_AMP_FILE'//itext, par_string=par%temp_amps(comp))
    do j = 1, par%numband
       call int2string(j,jtext)
       call get_parameter_hashtable(htbl, 'COMP'//trim(itext)//'_FIT'//jtext,&
            len_itext=len_jtext,par_lgt=par%fg_temp_corr(comp,j))
       if (par%fg_temp_corr(comp,j)) then
          par%fg_nfit(comp) = par%fg_nfit(comp) + 1
       end if
    end do
    call get_parameter_hashtable(htbl, 'COMP_POLTYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_poltype(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR_GAUSS_MEAN'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,1,1))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR_GAUSS_STD'//itext, len_itext=len_itext,&
         par_dp=par%fg_gauss(comp,1,2))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR_UNI_LOW'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,1,1))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR_UNI_HIGH'//itext, len_itext=len_itext,&
         par_dp=par%fg_uni(comp,1,2))
    call get_parameter_hashtable(htbl, 'COMP_T_POLTYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_poltype(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T'//itext, len_itext=len_itext, par_dp=par%fg_init(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_SAMP_NSIDE'//itext, len_itext=len_itext,&
         par_int=par%fg_samp_nside(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_SAMPLE'//itext, len_itext=len_itext,&
         par_lgt=par%fg_samp_spec(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_INPUT_MAP'//itext, len_itext=len_itext,&
         par_string=par%fg_spec_file(comp,1))       
    call get_parameter_hashtable(htbl, 'COMP_T_REGION'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_region(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_LNL_TYPE'//itext, len_itext=len_itext,&
         par_string=par%fg_ind_lnl(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_PRIOR'//itext, len_itext=len_itext,&
         par_string=par%fg_prior_type(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_STEPSIZE'//itext,len_itext=len_itext,&
         par_dp=par%fg_spec_step(comp,1))
    call get_parameter_hashtable(htbl, 'COMP_T_TUNE_STEPSIZE'//itext,len_itext=len_itext,&
         par_lgt=par%fg_spec_tune(comp,1))

  end subroutine read_hi_fit
    
  function get_token(string, sep, num, group, allow_empty) result(res)
    implicit none
    character(len=*)           :: string, sep
    character(len=len(string)) :: res
    character(len=*), optional :: group
    logical(lgt),     optional :: allow_empty
    integer(i4b)               :: i, num, ext(2)
    ext = -1
    do i = 1, num
       call tokenize(string, sep, ext, group, allow_empty)
    end do
    res = string(ext(1):ext(2))
  end function get_token
  
  ! Fill all tokens into toks, and the num filled into num                                           
  subroutine get_tokens(string, sep, toks, num, group, maxnum, allow_empty)
    implicit none
    character(len=*) :: string, sep
    character(len=*) :: toks(:)
    character(len=*), optional :: group
    integer(i4b),     optional :: num, maxnum
    logical(lgt),     optional :: allow_empty
    integer(i4b) :: n, ext(2), nmax
    ext = -1
    n = 0
    nmax = size(toks); if(present(maxnum)) nmax = maxnum
    call tokenize(string, sep, ext, group, allow_empty)
    do while(ext(1) > 0 .and. n < nmax)
       n = n+1
       toks(n) = string(ext(1):ext(2))
       call tokenize(string, sep, ext, group, allow_empty)
    end do
    if(present(num)) num = n
  end subroutine get_tokens
  
  subroutine tokenize(string, sep, ext, group, allow_empty)
    implicit none
    character(len=*) :: string, sep
    character(len=*), optional   :: group
    character(len=256)  :: op, cl
    integer(i4b), save           :: level(256), nl
    integer(i4b), intent(inout)  :: ext(2)
    logical(lgt), optional       :: allow_empty
    
    integer(i4b) :: i, j, o, c, ng
    logical(lgt) :: intok, hit, empty
    
    empty = .false.; if(present(allow_empty)) empty = allow_empty
    
    if(ext(2) >= len(string)) then
       ext = (/ 0, -1 /)
       return
    end if
    ng = 0
    if(present(group)) then
       ng = len_trim(group)/2
       do i = 1, ng
          op(i:i) = group(2*i-1:2*i-1)
          cl(i:i) = group(2*i:2*i)
       end do
    end if
    if(ext(2) <= 0) then
       level = 0
       nl = 0
    end if
    intok = .false.
    j     = 1
    do i = ext(2)+2, len(string)
       hit = .false.
       c = index(cl(1:ng), string(i:i))
       if(c /= 0) then; if(level(c) > 0) then
          level(c) = level(c) - 1
          if(level(c) == 0) nl = nl - 1
          hit = .true.
       end if; end if
       if(nl == 0) then
          ! Are we in a separator or not                                                             
          if(index(sep, string(i:i)) == 0) then
             ! Nope, so we must be in a token. Register start of token.                              
             if(.not. intok) then
                j = i
                intok = .true.
             end if
          else
             ! Yes. This either means that a token is done, and we should                            
             ! return it, or that we are waiting for a new token, in                                 
             ! which case do nothing.                                                                
             if(intok) then
                ext = (/ j, i-1 /)
                return
             elseif(empty) then
                ext = (/ i, i-1 /)
                return
             end if
          end if
       end if
       o = index(op(1:ng), string(i:i))
       if(o /= 0 .and. .not. hit) then
          if(level(o) == 0) nl = nl + 1
          level(o) = level(o) + 1
       end if
    end do
    ! Handle last token                                                                              
    if(intok) then
       ext = (/ j, i-1 /)
    elseif(empty) then
       ext = (/ i, i-1 /)
    else
       ext = (/ 0, -1 /)
    end if
  end subroutine tokenize
  
end module dang_param_mod
