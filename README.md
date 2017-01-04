# Allink
--------------- OVERVIEW -------------------

Collection of utilities based on two basics classes: 
Matematica and VarData.

Matematica) performs math operations on vectors and matrices 
for smoothing, interpolation, convolution, image processing...

VarData) manipulate a structure of points connected by links.

Addraw) openGL engine.

ElPoly) analyze mechanical properties of polymer and membrane 
like structures.

Addyn) perform molecular dynamics and Monte Carlo simulations 
and has a solver for 4th oder PDE. 

Avvis) perform all the operation of Matematica on different 
sets of data visualized on a Qt graphical interface.

DrImage) image manipulation on the Matematica filters.

The program is intended to use as less as possible external 
libs (optional: openGL, gsl, fftw, cgal, png, tiff, boost, 
MPI, Qt...).

--------------- INSTALL --------------------

Edit the file /include/MakeInclude.mk to set the flags

Type make
