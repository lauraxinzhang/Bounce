/**
 * \file    Constants.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    March, 2019
 *
 * \brief   Declares the structs for root finding.
 *
 * \NOTE    See header file for detailed explainations for each 
 *          function \n
 *
 *          Uses Numerical Recipe nr3.h data containers for Vector and Matrix
 *          functionalities; uses Vector.h for 3D vectors with arithmetic operations.
 *      
 *          Uses nr3 routines for interpolations.     
 *      
 *          Uses double instead of float to avoid floating number errors.
 */

#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED 1

// Physical contants. Don't change unless you find the value to be wrong.
#define PI 3.14159265358979323846		 // pi, no explanation needed...
#define CC 299792458.0                   // speed of light

#define ME 9.10938356E-31                       // electron mass, in kilogram
#define MI 1.6726219E-27                 // ion mass, in kilogram
#define QE 1.60217662E-19                // elementary charge, in coulomb
#define EPSILON0 8.8541878128E-12                // F/m SI units.

#define EVTOJOULE 1.60217662E-19         // energy conversion factor
#define TESLATOGAUSS 10000               // unit conversion

#define NPERORBIT 20
#define COULOMBLOG 19                    // coulomb logarithm for collisions TODO: calculate this self consistently
#define PRECISION 1E-15
#define INTERFACE 0.92                   // cutoff in psinorm that defines the core/edge interface

// Numbers specific to LTX geometry
//#define BMAGAXIS  2000                   // mod(B) = 0.2 T = 2000 Gauss at magnetic axis, characteristic field strength
                     // steps per Lamor orbit, from Boris convergence test.
//#define RMAJOR    0.4

// desired size of grid spacing
// Recommended values for conventional aspect ratio:
//#define DSMIN 0.1 // minimum size of triangles, in percentage of major radius
//#define DSMAX 2 // maximum size of triangles, in percentage of major radius


// Recommended values for ST:
#define DSMIN 2 // minimum size of triangles, in percentage of major radius
#define DSMAX 8 // maximum size of triangles, in percentage of major radius
#define COREZONE 50

#define MAXREFINESTEPS 5
#define RSEPOFFSET  0.00 // m

// which interpolation to use for magnetic field
// if USELINEAR = 1, linearly interpolated
// if USELINEAR = 0, boost library BSpline 
#define USELINEAR 1

#endif
