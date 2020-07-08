module param_mod
    use healpix_types
    use init_mod
    use hashtbl
    implicit none

    type params
        integer(i4b)   :: numband
        integer(i4b)   :: ncomp
        integer(i4b)   :: ngibbs
        character(len=512) :: datadir
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

        ! write(*,*) paramfile
        call get_file_length(paramfile,parfile_len)
        allocate(parfile_cache(parfile_len))
        call read_paramfile_to_ascii(paramfile,parfile_cache)

        !Initialize a hash table                                                                         
        call init_hash_tbl_sll(htable,tbl_len=10*parfile_len)
        ! Put the parameter file into the hash table                                                     
        call put_ascii_into_hashtable(parfile_cache,htable)
        deallocate(parfile_cache)
    
        

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
1           close(units(depth))
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
    
    ! subroutine read_global_params

    ! end subroutine read_global_params

    ! subroutine read_data_params

    ! end subroutine read_data_params

    ! subroutine read_comp_params

    ! end subroutine read_comp_params

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