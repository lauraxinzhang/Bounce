/**
 * \file    GEqlibrium.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    December, 2020
 *
 * \brief   Declares the GEquilibrium class. Parses the gEqdsk file, calculates and interpolates B field
 *         components and normalized poloidal flux
 */

#ifndef GE_H_INCLUDED
#define GE_H_INCLUDED 1

#include <cstdio>
#include <time.h> // for seeding the random engine.
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "Containers.hpp"
#include "Parser.hpp"
#include "Constants.hpp"
#include "Vector.hpp"
#include "Mesh.hpp"
#include "UniMesh.hpp"
#include "Fields.hpp"

class GEquilibrium : public Fields
{
public:
    /**
    *\brief Default constructor
    */
    GEquilibrium();
    
    /**
     * \brief Constructs data members according to given size
     * \note Used in constructing testing field configurations
     */
    GEquilibrium(int nw, int nh);
    
    /**
     *\brief Construct from geqdsk file
     *\param gpath path to eqdsk file
     */
    GEquilibrium(const std::string& gpath);
    
    /**
     * \brief Calculates the value of magnetic field on the input grid
     * \note Nothing. Values for B are stored in private members.
     */
    void setupB();
    
    /**
     *\brief Sets up the limiter and separatrix lineseg vectors
     */
    void setupLimiter();
    
    /**
     * \brief Interpolates on uniform grid and get B for given location p(r, z)
     * \return a Vector containing (Br, Btor, Bz)
     */
    Vector getB(MeshPoint& p);
    
    /**
     * \brief Interpolates on uniform grid and get B for given location v(x, y, z)
     * \param v Location of query in cartesian coordinates
     * \return a Vector containing (Bx, By, Bz)
     * \details the input position vector v is first converted to (r, z), then getB(MeshPoint& p) is called. The resulting
     * B vector in cylindrical geomatry is then converted to cartesian geometry and then returned.
     */
    Vector getB(Vector& v);
    
    /**
     * \brief Interpolates and finds psi at given location p(r, z)
     */
    double getPsi(MeshPoint& p);
    
    /**
     *\brief Get normalized psi at given location
     */
    double getPsiNorm(MeshPoint& p);
    
    /**
     *\brief Get normalized psi at given location
     */
    double getPsiNorm(Vector& v);
    
    /**
     *\brief get poloidal angle
     */
    double getCosTheta(Vector& v);
    
    /**
     * \brief Whether the location v is at the limiter
     * \note Must be inherited by child class
     */
    bool isLimiter(Vector v);
    
    /**
     * \brief Whether the location v is at the limiter
     * \param v Location to check
     * \param wallSeg Updated to the wall segment that the point intersects with
     * \note Must be inherited by child class
     */
    bool isLimiter(Vector v, LineSeg& wallSeg);
    
    /**
     * \brief Whether the location v is outside the separatrix
     */
    bool isSeparatrix(Vector v);
    
    /**
     * \brief Whether the location p is outside the separatrix
     */
    bool isSeparatrix(MeshPoint p);
    
    MatDoub Br_, Bz_, Btor_;
    VecDoub rGrid_, zGrid_, siGrid_;
    
    Parser::eqdsk gfile_;
    double rlimHFS_, rlimLFS_;
    double rsepHFS_, rsepLFS_;
    MeshPoint xtop_, xbot_; // xpoints
    MeshPoint maxis_; // magnetic axis
    
protected:
    UniMesh * mBr_;
    UniMesh * mBz_;
    UniMesh * mBtor_;
    UniMesh * mPsiRZ_;
    UniMesh * mRPsi_;
    std::vector<LineSeg> limiter_;
    std::vector<LineSeg> separatrix_;
    
    /**
     * \brief Coefficients for barycentric interpolation
     */
    int gamma(int k);
    
    /**
     * \brief Calculates differentiation matrix for rational interpolation
     * \param grid rGrid or zGrid, the orientation of interpolation
     * \note Interpolate along r if r derivative is needed, or z if z derivative is needed.
     */
    MatDoub Dij(VecDoub& grid);
    
    
    
    friend class FieldsTests;
    friend class Grid;
    
};


#endif
