c
c -----------------------------------------------------
c
      subroutine valout(lst, lend, time, nvar, naux)
c     subroutine valout_slice(lst, lend, time, nvar, naux)
c
c     # Output the results for a general system of conservation laws
c     # in 3 dimensions as a 2d amrclaw output file, along one 
c     # coordinate-aligned slice.

c     # As initial test, use slice z=zupper
c
c     # Write the results to the file fort.q<iframe>
c
c     # set outaux = .true. to also output the aux arrays to fort.a<iframe>

      use amr_module
c     use slice_module, only: icoord, xyz_slice

      implicit double precision (a-h,o-z)
      character*10  fname1, fname2, fname3, fname4, fname5

      logical outaux
      integer output_aux_num
      logical intersects
      double precision val(nvar)


      iadd(ivar,i,j,k)   = loc     +    (ivar-1)
     &                             +    (i-1)*nvar
     &                             +    (j-1)*nvar*mitot
     &                             +    (k-1)*nvar*mitot*mjtot
      iaddaux(iaux,i,j,k) = locaux +    (iaux-1)
     &                             +    (i-1)*naux
     &                             +    (j-1)*naux*mitot
     &                             +    (k-1)*naux*mitot*mjtot

c     # hard-wire z=zupper (from amr_module)
      icoord = 3  ! means z slice
c     xyz_slice = zupper - 1.d-3
      xyz_slice = 0.6d0
      tol = 1.d-6

      output_aux_num = 0
      do i = 1, naux
        output_aux_num = output_aux_num + output_aux_components(i)
      end do

c     # Currently outputs all aux components if any are requested!
      outaux = ((output_aux_num > 0) .and. 
     .         ((.not. output_aux_onlyonce) .or. (time==t0)))

c     open(unit=77,file='fort.b',status='unknown',access='stream')


c     ### Python graphics output
c
c        ###  make the file names and open output files
      fname1 = 'fort.qxxxx'
      fname2 = 'fort.txxxx'
      fname3 = 'fort.axxxx'
      fname4 = 'fort.bxxxx'
      matunit1 = 50
      matunit2 = 60
      matunit3 = 70
      matunit4 = 71
      nstp     = matlabu
      do 55 ipos = 10, 7, -1
         idigit = mod(nstp,10)
         fname1(ipos:ipos) = char(ichar('0') + idigit)
         fname2(ipos:ipos) = char(ichar('0') + idigit)
         fname3(ipos:ipos) = char(ichar('0') + idigit)
         fname4(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp / 10
 55   continue

      open(unit=matunit1,file=fname1,status='unknown',
     .       form='formatted')

      if (output_format == 3) then
c         # binary output          
          open(unit=matunit4,file=fname4,status='unknown',
     &            access='stream')
          endif

      intersects = .false.
      level = lst
      ngrids = 0
c65   if (level .gt. lfine) go to 90
 65   if (level .gt. lend) go to 90
         mptr = lstart(level)
 70      if (mptr .eq. 0) go to 80
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
           zup = zlow + nz*hzposs(level)
           write(6,*) '+++ xyz_slice, zlow, zup: '
           write(6,*) '+++ ',xyz_slice, zlow, zup
           if ((icoord==3) .and. (xyz_slice .ge. zlow - tol) .and.
     &             (xyz_slice .le. zup + tol)) then
              intersects = .true.
              ngrids  = ngrids + 1
              write(matunit1,1001) mptr, level, nx, ny
 1001         format(i5,'                 grid_number',/,
     &               i5,'                 AMR_level',/,
     &               i5,'                 mx',/,
     &               i5,'                 my')

              write(matunit1,1002) xlow,ylow,hxposs(level),
     &           hyposs(level)
 1002        format(e18.8,'    xlow', /,
     &              e18.8,'    ylow', /,
     &              e18.8,'    dx', /,
     &              e18.8,'    dy')
           endif


         write(6,*) '+++ mptr, intersects: ',mptr,intersects

         if ((output_format == 1) .and. intersects) then
            if (icoord==3) then
              do 75 k = nghost+1, mktot-nghost
                   zkm = zlow + (k-nghost-1)*hzposs(level)
                   zkp = zlow + (k-nghost)*hzposs(level)
                   alpha = -1.d0  !# indicates ignore this row
                   if ((xyz_slice .ge. zkm) .and. 
     &                 (xyz_slice .lt. zkp)) then
                      alpha = (xyz_slice - zkm)/hzposs(level)
                    else if ((k==nz) .and. (xyz_slice==zupper) .and.
     &                     (zkp > zupper - 0.5d0*hzposs(level))) then
c                     # special case for upper domain boundary
                      alpha = 1.d0
                    endif
                  
c                  write(6,*) '+++ k, zkm, alpha: ',k, zkm, alpha
                   if (alpha > -1.d0) then
                      do 76 j = nghost+1, mjtot-nghost
                      do 77 i = nghost+1, mitot-nghost
                      do ivar=1,nvar
                         val(ivar) = alpha*alloc(iadd(ivar,i,j,k+1)) 
     &                         + (1.d0-alpha)*alloc(iadd(ivar,i,j,k))
                         if (dabs(val(ivar)) < 1d-90) val(ivar) = 0.d0
                         enddo
                      write(matunit1,109) (val(ivar), ivar=1,nvar)

  109                 format(50e26.16)
   77                 continue
                      write(matunit1,*) ' '
   76                 continue
               endif
               write(matunit1,*) ' '
   75       continue
            write(matunit1,*) ' '
         endif
         endif

         if ((output_format == 3) .and. intersects) then
c          # binary output          
           write(6,*) 'binary output not implemented'
           endif


         mptr = node(levelptr, mptr)
         go to 70
   80 level = level + 1
      go to 65

   90 continue


c     --------------
c     # fort.t file:
c     --------------

      open(unit=matunit2,file=fname2,status='unknown',
     .       form='formatted')

      ndim = 3
c     # NOTE: we need to print out nghost too in order to strip
c     #       ghost cells from q when reading in pyclaw.io.binary
      write(matunit2,1000) time,nvar,ngrids,naux,2,nghost
 1000 format(e18.8,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 naux'/,
     &       i5,'                 ndim'/,
     &       i5,'                 nghost'/,/)
c

      write(6,601) matlabu,time
  601 format('AMRCLAW: Frame ',i4,
     &       ' output files done at time t = ', d12.6,/)

      matlabu = matlabu + 1

      close(unit=matunit1)
      close(unit=matunit2)
      if (output_format == 3) then
          close(unit=matunit4)
          endif

      return
      end
