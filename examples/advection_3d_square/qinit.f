c
c
c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,mz, xlower,ylower,zlower,
     &                  dx,dy,dz,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit double precision (a-h,o-z)
c
       dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
       dimension x(1-mbc:mx+mbc)
       dimension aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc)
c
c     # set concentration profile
c     ---------------------------
c

       do i = 1,mx
          xi = xlower + (i-0.5d0)*dx
             do j = 1,my
                yj = ylower + (j-0.5d0)*dy
                do k = 1,mz
                   zk = zlower + (k-0.5d0)*dz
                   if ((xi.gt.0.2d0) .and. (xi.lt.0.5d0) .and.
     &                 (yj.gt.0.2d0) .and. (yj.lt.0.5d0) .and. 
     &                 (zk.gt.0.2d0) .and. (zk.lt.0.5d0)) then
                          q(1,i,j,k) = 1.d0
                      else
                          q(1,i,j,k) = 0.d0
                      endif
                enddo
             enddo
         enddo

      return
      end
