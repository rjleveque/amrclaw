
.. _amrclaw_examples_advection_2d_square:

Two-dimensional advection of a square pulse 
===========================================

With periodic boundary conditions.

Modified to illustrate how to allow refinement only in some region even
after flagged cells are buffered.  

`flag2refine2.f90` calls `allowflag.f90` which checks this.
