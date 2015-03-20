module slices_module
!     # Output the results for a general system of conservation laws
!     # in 3 dimensions as a 2d amrclaw output file, along one 
!     # coordinate-aligned slice.


    implicit none
    save
    
    ! Number of slices and their position
    integer :: num_slices, num_slices_xy, num_slices_xz, num_slices_yz
    real (kind=8), allocatable, dimension(:) :: slices_xy, slices_xz, slices_yz
    
    ! Tolerance for deciding whether point is in slice
    real (kind=8), parameter :: tol = 1.d-6
    
    
contains
    
    subroutine set_slices(fname)
    
        implicit none
        
        ! Input
        character(len=*), intent(in), optional :: fname
        
        ! Local variables
        integer :: i
        integer :: iunit = 7
        
        ! Open file
        if (present(fname)) then
            call opendatafile(iunit, fname)
        else
            call opendatafile(iunit, 'slices.data')
        endif
        
        read(iunit,*) num_slices_xy
        read(iunit,*) num_slices_xz
        read(iunit,*) num_slices_yz

        num_slices = num_slices_xy + num_slices_xz + num_slices_yz
        
        allocate(slices_xy(num_slices_xy), slices_xz(num_slices_xz), slices_yz(num_slices_yz))
        
        ! Read in xy-slice position in order
        do i=1,num_slices_xy
            read(iunit,*) slices_xy(i)
        end do
        
        ! Read in xz-slice position in order
        do i=1,num_slices_xz
            read(iunit,*) slices_xz(i)
        end do
        
        ! Read in yz-slice position in order
        do i=1,num_slices_yz
            read(iunit,*) slices_yz(i)
        end do
        
        write(6,*) slices_xy(1)
        
        close(iunit)
        
    end subroutine set_slices

! ########################################################################
    
    subroutine print_slices_xy(lst, lend, time, nvar, naux)

        use amr_module

        integer :: lst, lend, naux, nvar
        real (kind=8) :: time
        
        character*15 :: fname1, fname2
        real (kind=8) :: val(nvar)
        real (kind=8) :: alpha, xlow, xup, ylow, yup, zlow, zup, zkm, zkp
        
        integer :: i, j, k, l, ipos, ivar, iaux, idigit, kprint, level, loc, locaux
        integer :: matunit1, matunit2, mitot, mjtot, mktot, mptr, ndim, ngrids
        integer :: nx, ny, nz, nstp
        integer iadd
        integer iaddaux
        logical in_grid
        
        iadd(ivar,i,j,k)   = loc     +    (ivar-1) &
                                  +    (i-1)*nvar &
                                  +    (j-1)*nvar*mitot &
                                  +    (k-1)*nvar*mitot*mjtot
        iaddaux(iaux,i,j,k) = locaux +    (iaux-1) &
                                  +    (i-1)*naux &
                                  +    (j-1)*naux*mitot &
                                  +    (k-1)*naux*mitot*mjtot
        
                
!     ###  make the file names and open output files
        fname1 = 'slice_xyx.qxxxx'
        fname2 = 'slice_xyx.txxxx'
!        fname3 = 'fort.axxxx'
!        fname4 = 'fort.bxxxx'
        matunit1 = 50
        matunit2 = 60
 !       matunit3 = 70
 !       matunit4 = 71
        
        nstp     = matlabu
        do ipos = 15, 12, -1
            idigit = mod(nstp,10)
            fname1(ipos:ipos) = char(ichar('0') + idigit)
            fname2(ipos:ipos) = char(ichar('0') + idigit)
!            fname3(ipos:ipos) = char(ichar('0') + idigit)
!            fname4(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
        end do
        
        do i = 1,num_slices_xy
            fname1(9:9) = char(ichar('0') + i)
            open(unit=matunit1+i,file=fname1,status='unknown',form='formatted')
        end do

!      if (output_format == 3) then
!         # binary output          
!          open(unit=matunit4,file=fname4,status='unknown',
!     &            access='stream')
!          endif

        level = lst
        ngrids = 0

        do while (level .le. lend)
            mptr = lstart(level)
                    
            do while (mptr .ne. 0)
                nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
                ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
                nz      = node(ndkhi,mptr) - node(ndklo,mptr) + 1
                loc     = node(store1, mptr)
                locaux  = node(storeaux,mptr)
                mitot   = nx + 2*nghost
                mjtot   = ny + 2*nghost
                mktot   = nz + 2*nghost
                xlow = rnode(cornxlo,mptr)
                ylow = rnode(cornylo,mptr)
                zlow = rnode(cornzlo,mptr)
                xup = xlow + nx*hxposs(level)
                yup = ylow + ny*hyposs(level)
                zup = zlow + nz*hzposs(level)

                do l = 1,num_slices_xy

                    in_grid = .false.   ! initialize for this grid and check below

                    if ((slices_xy(l) .ge. zlow - tol) .and. (slices_xy(l) .le. zup + tol)) then
                        in_grid = .true.
                        ngrids  = ngrids + 1
                        write(matunit1+l,1301) mptr, level, nx, ny
           1301         format(i5,'                 grid_number',/, &
                               i5,'                 AMR_level',/, &
                               i5,'                 mx',/, &
                               i5,'                 my')
                        
                        write(matunit1+l,1302) xlow, ylow, hxposs(level), hyposs(level)
            1302        format(e18.8,'    xlow', /, &
                               e18.8,'    ylow', /, &
                               e18.8,'    dx', /, &
                               e18.8,'    dy')                           
                    end if

                    if ((output_format == 1) .and. in_grid) then
                        write(matunit1+l,*) ' '
                        kprint = -1
                        do k = nghost+1, mktot-nghost
                            zkm = zlow + (k-nghost-1)*hzposs(level)
                            zkp = zlow + (k-nghost)*hzposs(level)
                            alpha = -1.d0  !# indicates ignore this row
                            if ((slices_xy(l) .ge. zkm) .and. (slices_xy(l) .lt. zkp)) then
                                alpha = (slices_xy(l) - zkm)/hzposs(level)
                            else if ((k==nghost+1) .and. (slices_xy(l) .le. zkm)) then
    !                     # special case for lower edge
                                alpha = 0.d0
                            else if ((k==mktot-nghost) .and. (slices_xy(l)==zupper) .and. &
                                    (zkp > zupper - 0.5d0*hzposs(level))) then
    !                     # special case for upper domain boundary
                                alpha = 1.d0
                            endif
                      
                            if (alpha > -1.d0) then
                                kprint = k
                                do j = nghost+1, mjtot-nghost
                                    do i = nghost+1, mitot-nghost
                                        do ivar=1,nvar
                                            val(ivar) = alpha*alloc(iadd(ivar,i,j,k+1)) & 
                                                 + (1.d0-alpha)*alloc(iadd(ivar,i,j,k))
                                            if (dabs(val(ivar)) < 1d-90) then 
                                                val(ivar) = 0.d0
                                            end if
                                        end do
                                        write(matunit1+l,109) (val(ivar), ivar=1,nvar)
                    109                 format(50e26.16)
                                    end do ! i loop
                                    write(matunit1+l,*) ' '
                                end do ! j loop
                            endif
                            write(matunit1+l,*) ' '                            
                        end do ! k loop
                    
                        if (kprint < 0) then
                            write(6,*) '*** ERROR in slices_module'
                            write(6,*) 'slices_xy: ', slices_xy(l)
                            write(6,*) 'zlow, zup, dz: ',zlow,zup,hzposs(level)
                            stop
                        end if
                        write(matunit1+l,*) ' '
                    end if
                end do
            
!         if ((output_format == 3) .and. intersects) then
!          # binary output          
!          write(6,*) 'binary output not implemented'
!          endif

                mptr = node(levelptr, mptr)                
            end do
            level = level + 1
        end do

    !     --------------
    !     # fort.t file:
    ! --------------




        do i = 1,num_slices_xy
            fname2(9:9) = char(ichar('0') + i)
            open(unit=matunit2+i,file=fname2,status='unknown',form='formatted')
            ndim = 3
        !     # NOTE: we need to print out nghost too in order to strip
        !     #       ghost cells from q when reading in pyclaw.io.binary
            write(matunit2+i,1000) time,nvar,ngrids,naux,2,nghost,slices_xy(i)
            1000 format(e18.8,'    time', /, &
                i5,'                 meqn'/, &
                i5,'                 ngrids'/, &
                i5,'                 naux'/, &
                i5,'                 ndim'/, &
                i5,'                 nghost'/, &
                e18.8,'    location'/)
            
            close(unit=matunit1+i)
            close(unit=matunit2+i)
            !      if (output_format == 3) then
            !         close(unit=matunit4)
            !       endif
        end do
        
    end subroutine print_slices_xy 

    
! #######################################################################

    subroutine print_slices_yz(lst, lend, time, nvar, naux)

        use amr_module

        integer :: lst, lend, naux, nvar
        real (kind=8) :: time
        
        character*15 :: fname1, fname2, fname3, fname4
        real (kind=8) :: val(nvar)
        real (kind=8) :: alpha, xlow, xup, ylow, yup, zlow, zup, xim, xip
        
        integer :: i, j, k, l, ipos, ivar, iaux, idigit, iprint, level, loc, locaux
        integer :: matunit1, matunit2, mitot, mjtot, mktot, mptr, ndim, ngrids
        integer :: nx, ny, nz, nstp
        integer iadd
        integer iaddaux
        logical in_grid
        
        iadd(ivar,i,j,k)   = loc     +    (ivar-1) &
                                  +    (i-1)*nvar &
                                  +    (j-1)*nvar*mitot &
                                  +    (k-1)*nvar*mitot*mjtot
        iaddaux(iaux,i,j,k) = locaux +    (iaux-1) &
                                  +    (i-1)*naux &
                                  +    (j-1)*naux*mitot &
                                  +    (k-1)*naux*mitot*mjtot
        

!     ###  make the file names and open output files
        fname1 = 'slice_yzx.qxxxx'
        fname2 = 'slice_yzx.txxxx'
!        fname3 = 'fort.axxxx'
!        fname4 = 'fort.bxxxx'
        matunit1 = 50
        matunit2 = 60
 !       matunit3 = 70
 !       matunit4 = 71
        
        nstp     = matlabu
        do ipos = 15, 12, -1
            idigit = mod(nstp,10)
            fname1(ipos:ipos) = char(ichar('0') + idigit)
            fname2(ipos:ipos) = char(ichar('0') + idigit)
!            fname3(ipos:ipos) = char(ichar('0') + idigit)
!            fname4(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
        end do
        
        do i = 1,num_slices_yz
            fname1(9:9) = char(ichar('0') + i)
            open(unit=matunit1+i,file=fname1,status='unknown',form='formatted')
        end do

!      if (output_format == 3) then
!         # binary output          
!          open(unit=matunit4,file=fname4,status='unknown',
!     &            access='stream')
!          endif

        level = lst
        ngrids = 0

        do while (level .le. lend)
            mptr = lstart(level)
                    
            do while (mptr .ne. 0)
                nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
                ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
                nz      = node(ndkhi,mptr) - node(ndklo,mptr) + 1
                loc     = node(store1, mptr)
                locaux  = node(storeaux,mptr)
                mitot   = nx + 2*nghost
                mjtot   = ny + 2*nghost
                mktot   = nz + 2*nghost
                xlow = rnode(cornxlo,mptr)
                ylow = rnode(cornylo,mptr)
                zlow = rnode(cornzlo,mptr)
                xup = xlow + nx*hxposs(level)
                yup = ylow + ny*hyposs(level)
                zup = zlow + nz*hzposs(level)


                do l = 1,num_slices_yz

                    in_grid = .false.   ! initialize for this grid and check below
    
                    if ((slices_yz(l) .ge. xlow - tol) .and. (slices_yz(l) .le. xup + tol)) then
                    
                        in_grid = .true.
                        ngrids  = ngrids + 1
                        write(matunit1+l,1301) mptr, level, ny, nz
           1301         format(i5,'                 grid_number',/, &
                               i5,'                 AMR_level',/, &
                               i5,'                 my',/, &
                               i5,'                 mz')
                        
                        write(matunit1+l,1302) ylow, zlow, hyposs(level), hzposs(level)
            1302        format(e18.8,'    ylow', /, &
                               e18.8,'    zlow', /, &
                               e18.8,'    dy', /, &
                               e18.8,'    dz')                           
                    end if

                    if ((output_format == 1) .and. in_grid) then
                        write(matunit1+l,*) ' '
                        iprint = -1
                        do i = nghost+1, mitot-nghost
                            xim = xlow + (i-nghost-1)*hxposs(level)
                            xip = xlow + (i-nghost)*hxposs(level)
                            alpha = -1.d0  !# indicates ignore this row
                            if ((slices_yz(l) .ge. xim) .and. (slices_yz(l) .lt. xip)) then
                                alpha = (slices_yz(l) - xim)/hxposs(level)
                            else if ((i==nghost+1) .and. (slices_yz(l) .le. xim)) then
    !                     # special case for lower edge
                                alpha = 0.d0
                            else if ((i==mitot-nghost) .and. (slices_yz(l) .ge. xip)) then
    !                     # special case for upper domain boundary
                                alpha = 1.d0
                            endif
                      
                            if (alpha > -1.d0) then
                                iprint = i
                                do k = nghost+1, mktot-nghost
                                    do j = nghost+1, mjtot-nghost
                                        do ivar=1,nvar
                                            val(ivar) = alpha*alloc(iadd(ivar,i+1,j,k)) & 
                                                 + (1.d0-alpha)*alloc(iadd(ivar,i,j,k))
                                            if (dabs(val(ivar)) < 1.d-90) then 
                                                val(ivar) = 0.d0
                                            end if
                                        end do
                                        write(matunit1+l,109) (val(ivar), ivar=1,nvar)
                    109                 format(50e26.16)
                                    end do ! j loop
                                    write(matunit1+l,*) ' '
                                end do ! k loop
                            endif
                        end do ! i loop
                    
                        if (iprint < 0) then
                            write(6,*) '*** ERROR in slices_module'
                            write(6,*) 'slices_yz: ', slices_yz(l)
                            write(6,*) 'xlow, xup, dx: ',xlow,xup,hxposs(level)
                            stop
                        end if
                        write(matunit1+l,*) ' '
                    end if
                end do
            
!         if ((output_format == 3) .and. intersects) then
!          # binary output          
!          write(6,*) 'binary output not implemented'
!          endif

                mptr = node(levelptr, mptr)
            end do
            level = level + 1
        end do

    !     --------------
    !     # fort.t file:
    ! --------------

        do i = 1,num_slices_yz
            fname2(9:9) = char(ichar('0') + i)
            open(unit=matunit2+i,file=fname2,status='unknown',form='formatted')
            ndim = 3
        !     # NOTE: we need to print out nghost too in order to strip
        !     #       ghost cells from q when reading in pyclaw.io.binary
            write(matunit2+i,1000) time,nvar,ngrids,naux,2,nghost,slices_yz(i)
            1000 format(e18.8,'    time', /, &
                i5,'                 meqn'/, &
                i5,'                 ngrids'/, &
                i5,'                 naux'/, &
                i5,'                 ndim'/, &
                i5,'                 nghost'/, &
                e18.8,'    location'/)
                        
            close(unit=matunit1+i)
            close(unit=matunit2+i)
            !      if (output_format == 3) then
            !         close(unit=matunit4)
            !       endif
        end do
        
    end subroutine print_slices_yz

    
! #######################################################################


    subroutine print_slices(lst, lend, time, nvar, naux)
    
        integer :: lst, lend, naux, nvar
        real (kind=8) :: time

        if (num_slices_xy .gt. 0) then
            call print_slices_xy(lst, lend, time, nvar, naux)
        end if
!        if (num_slices_xz .gt. 0) print_slices_xz(lst, lend, time, nvar, naux)
        if (num_slices_yz .gt. 0) then 
            call print_slices_yz(lst, lend, time, nvar, naux)
        end if

    end subroutine print_slices
        
        
        
    
            
    
    
end module slices_module
    
    
