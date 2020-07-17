# reconstruction.jl
This repository contains code for the computer algebra system OSCAR developed within the TRR 195. 
It contains a very early version of Oscar code related to my PhD thesis at the University of Kaiserslautern.
Given an ideal $I$ in a polynomial ring over the rational numbers, it uses a Las Vegas algorithm to check whether
the given algebra defined by the ideal is isomorphic to the Tjurina algebra of a quasi-homogeneous isolated
hypersurface singularity or not. 

In case the result is positive it returns the defining equation $g$ of the hypersurface singularity and an automorphism $f$,
such that $f(I)$ is the Jacobian ideal of $g$.

For the moment you need a LINUX operating system and Julia 1.4.x or higher with a develpoment-version of 
OSCAR v0.4 to make it work.

How to install OSCAR and the code:
- Install JULIA 1.x.y with x => 4 on your computer (see https://github.com/JuliaLang/julia/tree/v1.4.2)
- Start JULIA
- Press "]" 
- run the command: 
```
dev "https://github.com/oscar-system/Oscar.jl" 
```
- Copy the folder "Singularities" and the file "reconstruction.jl" into "$HOME/.julia/dev/Oscar/examples/"

Here is an example:

After starting JULIA run the commands:
```
using Oscar
Oscar.example("reconstruction.jl")
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
I = ideal([-24*x^2*z+24*x*y*z-6*y^2*z+5*x^4, 12*x^2*z-12*x*y*z+3*y^2*z, -8*x^3+12*x^2*y-6*x*y^2+y^3+4*z^3]);
reconstruction(I)
J = ideal([4*x^3+x^2*y-5*x*y^2+2*y^3-3*x^2*z-3*x*y*z+4*y^2*z+2*x*z^2+3*y*z^2-2*z^3, x^3+x^2*y-5*x*y^2-5*y^3-4*x^2*z+5*x*y*z-3*y^2*z+5*x*z^2+y*z^2+3*z^3, -5*x^3-5*x*y^2+3*y^3-2*x*y*z-2*y^2*z+4*x*z^2+2*y*z^2-2*z^3]);
reconstruction(J)
```
