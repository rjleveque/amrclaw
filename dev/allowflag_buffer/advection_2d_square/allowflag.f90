! ::::::::::::::::::::: allowflag ::::::::::::::::::::::::::::::::::
!
!  Indicate whether the grid point at (x,y,t) at this refinement level
!  is allowed to be flagged for further refinement.
!
!  This is useful if you wish to zoom in on some structure in a 
!  known location but don't want the same level of refinement elsewhere.  
!  Points are flagged only if one of the errors is greater than the 
!  corresponding tolerance.
!
!  For example, to allow refinement of Level 1 grids everywhere but
!  of finer grids only for  y >= 0.4:
!  allowed(x,y,t,level) = (level.le.1 .or. y.ge.0.4d0) 
!
!  This routine is called from routine flag2refine.
!  If Richardson error estimates are used (if flag_richardson is true) 
!  then this routine is also called from errf1.


logical function allowflag(x,y,t,level)

    use amr_module, only: ibuff, hxposs, hyposs
    implicit none

    real(kind=8), intent(in) :: x,y,t
    integer, intent(in) :: level

    real (kind=8) :: x1,x2,y1,y2

    ! default version allows refinement anywhere:
    !allowflag = .true.

    x1 = 0.2d0
    x2 = 0.7d0
    y1 = 0.3d0
    y2 = 0.8d0

    allowflag = (x > x1 + ibuff*hxposs(level)) .and. &
                (x < x2 - ibuff*hxposs(level)) .and. &
                (y > y1 + ibuff*hyposs(level)) .and. &
                (y < y2 - ibuff*hyposs(level))

    end function allowflag
