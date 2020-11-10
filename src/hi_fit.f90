
    ! This architecture of this function has not been verified yet
    subroutine sample_HI_T(self, map_n) !(self, band, nside2, map_n)
        implicit none
  
        class(params)                                          :: self
        integer(i4b), intent(in)                               :: map_n!, nside2
        integer(i4b)                                           :: nside1, npix2, nside2
        ! real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: band
        real(dp), dimension(0:npix-1,nmaps,nbands)             :: cov
        real(dp), dimension(0:npix-1,nmaps)                    :: te
        real(dp), dimension(0:npix-1)                          :: te_sample
        real(dp), allocatable, dimension(:,:,:)                :: maps_low, cov_low
        real(dp), allocatable, dimension(:,:)                  :: T_low
        real(dp), allocatable, dimension(:)                    :: sample_T_low
        real(dp), dimension(nbands)                            :: signal, tmp
        real(dp), dimension(2)                                 :: x
        real(dp)                                               :: a, b, c, num, sam, t, p, sol, naccept_t_d

        te      = T_d
        cov     = dang_data%rms_map*dang_data%rms_map

        nside2 = 8

        nside1 = npix2nside(npix)
        if (nside1 == nside2) then
            npix2 = npix
        else
            npix2 = nside2npix(nside2)
        end if
        allocate(maps_low(0:npix2-1,nmaps,nbands))
        allocate(T_low(0:npix2-1,nmaps),cov_low(0:npix2-1,nmaps,nbands))
        allocate(sample_T_low(0:npix2-1))

        ! if (nside1 /= nside2) then 
        !     if (ordering == 1) then
        !         call udgrade_ring(te,nside1,T_low,nside2)
        !     else
        !         call udgrade_nest(te,nside1,T_low,nside2)
        !     end if
        !     do j = 1, nbands
        !         if (ordering == 1) then
        !             call udgrade_ring(map2fit(:,:,j),nside1,band_low(:,:,j),nside2)
        !             call convert_nest2ring(nside1,map2fit(:,:,j))
        !             call udgrade_ring(cov(:,:,j),nside1,cov_low(:,:,j),nside2)
        !             call convert_nest2ring(nside1,dang_data%rms_map(:,:,j))
        !         else
        !             call udgrade_nest(map2fit(:,:,j),nside1,band_low(:,:,j),nside2)
        !             call udgrade_nest(dang_data%rms_map(:,:,j),nside1,cov_low(:,:,j),nside2)
        !         end if
        !     end do
        !     cov_low = sqrt(cov_low / (npix/npix2))
        ! else
        do j = 1, nbands
            maps_low(:,:,j)   = dang_data%sig_map(:,:,j)
            cov_low(:,:,j)    = cov(:,:,j)
        end do
        T_low = te
        ! end if

        x(1) = 1.d0
        do i = 0, npix2-1
            if (mask(i,1) == 0.d0  .or. mask(i,1) == missval) then
                sample_T_low(i) = missval
                cycle
            else
                a   = 0.d0
                sol = T_low(i,map_n)

                ! Chi-square from the most recent Gibbs chain update
                do j = 1, nbands
                    a = a + (((temp01_amps(j,map_n)* HI(i,1)*planck(par%dat_nu(j)*1.d9,sol)) &
                        - maps_low(i,map_n,j))**2.d0)/cov_low(i,map_n,j)
                end do
                c   = a

                do l = 1, iterations
                    ! Begin sampling from the prior
                    ! t = rand_normal(sol,par%HI_Td_std)
                    t = rand_normal(par%HI_Td_mean,par%HI_Td_std)
                    b = 0.d0
                    do j = 1, nbands
                        tmp(j) = temp01_amps(j,map_n)*HI(i,1)*planck(par%dat_nu(j)*1.d9,t)
                        b      = b + ((tmp(j)-maps_low(i,map_n,j))**2.d0)/cov_low(i,map_n,j)
                    end do
                    b = b

                    if (b < c .and. t .lt. 35.d0 .and. t .gt. 10.d0) then
                        sam = t
                        c   = b
                    else
                        x(2) = exp(0.5d0*(c-b))
                        p = minval(x)
                        call RANDOM_NUMBER(num)
                        if (num < p) then
                            if (t .lt. 35.d0 .and. t .gt. 10.d0) then
                                sam = t
                                c   = b
                            end if
                        end if
                    end if
                end do
                sol             = sam
                sample_T_low(i) = sol
            end if
        end do
        if (nside1 /= nside2) then
            if (ordering == 1) then
                call udgrade_ring(sample_T_low, nside2, te_sample, nside1)
                call convert_nest2ring(nside2, sample_T_low)
            else
                call udgrade_nest(sample_T_low, nside2, te_sample, nside1)
            end if
        else
            te_sample =  sample_T_low
        end if
        T_d(:,1) = te_sample

        deallocate(maps_low)
        deallocate(T_low)
        deallocate(cov_low)
        deallocate(sample_T_low)

    end subroutine sample_HI_T
    ! ------------------------------------------------------------
