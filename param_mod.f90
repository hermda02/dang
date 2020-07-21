module param_mod
    use healpix_types
    use init_mod
    use utility_mod
    use hashtbl
    implicit none

    type params

        ! Global parameters
        integer(i4b)       :: ngibbs
        integer(i4b)       :: nsample  ! For things like the metrop-hast alg
        integer(i4b)       :: iter_out ! Out put maps every <- iterations
        logical(lgt)       :: output_fg
        character(len=512) :: outdir

        ! Data parameters
        integer(i4b)   :: numband
        character(len=512) :: datadir
        character(len=512), allocatable, dimension(:)   :: dat_label
        character(len=512), allocatable, dimension(:)   :: dat_mapfile
        character(len=512), allocatable, dimension(:)   :: dat_noisefile
        real(dp),           allocatable, dimension(:)   :: dat_nu

        ! Component parameters
        integer(i4b)   :: ncomp                                             ! # of foregrounds
        logical(lgt),       allocatable, dimension(:)     :: fg_inc         ! Logical - include fg?
        logical(lgt),       allocatable, dimension(:,:)   :: fg_sample_spec ! Logical - sample spec params
        logical(lgt),       allocatable, dimension(:)     :: fg_sample_amp  ! Logical - sample amplitude
        character(len=512), allocatable, dimension(:)     :: fg_label       ! Fg label (for outputs)
        character(len=512), allocatable, dimension(:)     :: fg_type        ! Fg type (power-law feks)
        real(dp),           allocatable, dimension(:)     :: fg_nu_ref      ! Fg reference frequency
        integer(i4b),       allocatable, dimension(:)     :: fg_ref_loc     ! Fg reference band
        real(dp),           allocatable, dimension(:,:,:) :: fg_gauss       ! Fg gaussian sampling
        real(dp),           allocatable, dimension(:,:,:) :: fg_uni         ! Fg sampling bounds
    end type params

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
            !write(*,*) "Exiting file " // filenames(depth)                                               
            depth = depth-1
        end do
        return

    end subroutine get_file_length

    subroutine read_param_file(par)
        implicit none
        type(hash_tbl_sll)                            :: htable
        type(params), intent(inout)                   :: par
        integer(i4b)                                  :: parfile_len, i
        character(len=512)                            :: paramfile
        character(len=512), allocatable, dimension(:) :: parfile_cache

        call getarg(1,paramfile)

        call get_file_length(paramfile,parfile_len)
        allocate(parfile_cache(parfile_len))
        call read_paramfile_to_ascii(paramfile,parfile_cache)

        !Initialize a hash table                                                                         
        call init_hash_tbl_sll(htable,tbl_len=10*parfile_len)
        ! Put the parameter file into the hash table                                                     
        call put_ascii_into_hashtable(parfile_cache,htable)
        deallocate(parfile_cache)
 
        call read_global_params(htable,par)    
        call read_data_params(htable,par)
        call read_comp_params(htable,par)

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
        ! write(*,*) units(depth)
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
        1   close(units(depth))
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
           if (key=="") cycle !we don't need blank lines                                                                                      
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

        ! logical(lgt)               :: found

        ! found = .false.
        ! call get_parameter_arg(parname, par_int, par_char, par_string, par_sp, par_dp, par_lgt, found, desc)
        ! if(found) then
        !     if(present(par_present)) par_present = .true.
        ! else
        call get_parameter_from_hash(htbl, parname, len_itext, par_int, &
            & par_char, par_string, par_sp, par_dp, par_lgt, par_present, desc)
        ! end if
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

        key=trim(parname)
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
        
        !if (cpar%myid == cpar%root) then                                                                
         
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
        type(params),       intent(inout) :: par

        integer(i4b)     :: i, j, n, len_itext
        character(len=2) :: itext
        character(len=2) :: jtext

        call get_parameter_hashtable(htbl, 'OUTPUT_DIRECTORY', par_string=par%outdir)
        call get_parameter_hashtable(htbl, 'NUMGIBBS', par_int=par%ngibbs)
        call get_parameter_hashtable(htbl, 'NUMSAMPLE', par_int=par%nsample)
        call get_parameter_hashtable(htbl, 'OUTPUT_ITER', par_int=par%iter_out)
        call get_parameter_hashtable(htbl, 'OUTPUT_COMPS', par_lgt=par%output_fg)
        
        par%outdir = trim(par%outdir) // '/'

    end subroutine read_global_params

    subroutine read_data_params(htbl,par)
        implicit none

        type(hash_tbl_sll), intent(in)    :: htbl
        type(params),       intent(inout) :: par

        integer(i4b)     :: i, j, n, len_itext
        character(len=2) :: itext
        character(len=2) :: jtext

        len_itext = len(trim(itext))

        call get_parameter_hashtable(htbl, 'NUMBAND',    par_int=par%numband)
        call get_parameter_hashtable(htbl, 'DATA_DIRECTORY', par_string=par%datadir)
        
        n = par%numband

        allocate(par%dat_mapfile(n),par%dat_label(n))
        allocate(par%dat_noisefile(n),par%dat_nu(n))

        do i = 1, n
            call int2string(i, itext)
            call get_parameter_hashtable(htbl, 'BAND_LABEL'//itext, len_itext=len_itext, par_string=par%dat_label(i))
            call get_parameter_hashtable(htbl, 'BAND_FILE'//itext, len_itext=len_itext, par_string=par%dat_mapfile(i))
            call get_parameter_hashtable(htbl, 'BAND_RMS'//itext, len_itext=len_itext, par_string=par%dat_noisefile(i))
            call get_parameter_hashtable(htbl, 'BAND_FREQ'//itext, len_itext=len_itext, par_dp=par%dat_nu(i))
        end do

    end subroutine read_data_params

    subroutine read_comp_params(htbl,par)
        implicit none

        type(hash_tbl_sll), intent(in)    :: htbl
        type(params),       intent(inout) :: par

        integer(i4b)     :: i, j, n, len_itext
        character(len=2) :: itext
        character(len=2) :: jtext

        len_itext = len(trim(itext))

        call get_parameter_hashtable(htbl, 'NUMCOMPS', par_int=par%ncomp)
        n = par%ncomp

        allocate(par%fg_label(n),par%fg_type(n),par%fg_nu_ref(n),par%fg_ref_loc(n))
        allocate(par%fg_inc(n),par%fg_sample_spec(n,2),par%fg_sample_amp(n))
        allocate(par%fg_gauss(n,2,2),par%fg_uni(n,2,2))

        do i = 1, n
            call int2string(i, itext)
            call get_parameter_hashtable(htbl, 'COMP_LABEL'//itext, len_itext=len_itext, par_string=par%fg_label(i))
            call get_parameter_hashtable(htbl, 'COMP_TYPE'//itext, len_itext=len_itext, par_string=par%fg_type(i))
            call get_parameter_hashtable(htbl, 'COMP_REF_FREQ'//itext, len_itext=len_itext, par_dp=par%fg_nu_ref(i))
            call get_parameter_hashtable(htbl, 'COMP_INCLUDE'//itext, len_itext=len_itext, par_lgt=par%fg_inc(i))

            if (trim(par%fg_type(i)) == 'power-law') then
               call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_BETA_MEAN'//itext, len_itext=len_itext,&
                    par_dp=par%fg_gauss(i,1,1))
               call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_BETA_STD'//itext, len_itext=len_itext,&
                    par_dp=par%fg_gauss(i,1,2))
               call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_BETA_LOW'//itext, len_itext=len_itext,&
                    par_dp=par%fg_uni(i,1,1))
               call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_BETA_HIGH'//itext, len_itext=len_itext,&
                    par_dp=par%fg_uni(i,1,2))
            else if (trim(par%fg_type(i)) == 'mbb') then
               call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_BETA_MEAN'//itext, len_itext=len_itext,&
                    par_dp=par%fg_gauss(i,1,1))
               call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_BETA_STD'//itext, len_itext=len_itext,&
                    par_dp=par%fg_gauss(i,1,2))
               call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_BETA_LOW'//itext, len_itext=len_itext,&
                    par_dp=par%fg_uni(i,1,1))
               call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_BETA_HIGH'//itext, len_itext=len_itext,&
                    par_dp=par%fg_uni(i,1,2))
               call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_T_MEAN'//itext, len_itext=len_itext,&
                    par_dp=par%fg_gauss(i,2,1))
               call get_parameter_hashtable(htbl, 'COMP_PRIOR_GAUSS_T_STD'//itext, len_itext=len_itext,&
                    par_dp=par%fg_gauss(i,2,2))
               call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_T_LOW'//itext, len_itext=len_itext,&
                    par_dp=par%fg_uni(i,2,1))
               call get_parameter_hashtable(htbl, 'COMP_PRIOR_UNI_T_HIGH'//itext, len_itext=len_itext,&
                    par_dp=par%fg_uni(i,2,2))
            end if
            par%fg_ref_loc(i) = minloc(abs(par%dat_nu-par%fg_nu_ref(i)),1)
        end do

    end subroutine read_comp_params

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

end module param_mod
