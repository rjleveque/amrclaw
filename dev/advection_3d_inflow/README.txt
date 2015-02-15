
Testing $CLAW/amrclaw/src/3d/valout_slice.f, which gives 2d-formatted
output on one coordinate slice rather than full 3d output.

Notes:

Makefile points to this version of valout.

Currently valout_slice.f has hardwired the coordinate direction and slice
location, and only has options icoord==1 (x-slice) or icoord==3 (z-slice).

E.g. if icoord==1 then xyz_slice is the x-location of the slice, which is in
the y-z plane.

Only one plane can be output currently in each run.

The setplot.py is set up as it would be for a 2d problem.


Future work:

Allow the user to specify a collection of planes in setrun.py, with values
of icoord and xyz_slice for each, and modify valout_slice.f to create a
separate output directory for each (Is this the best way? Or different
filenames for each in a single outdir?)

Allow specifying planes not aligned with coordinate axes?  This would be
more work to interpolate.  Note this is also needed if the user wants a
coordinate plane for a problem with a mapped grid.  Currently on a mapped
grid, xyz_slice would refer to a fixed value in the computational domain,
not the physical domain.


