/**
 * \file    Fields.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2020
 *
 * \brief   Declares the Fields class. Parses the gEqdsk file, calculates and interpolates B field
 *         components and normalized poloidal flux
 */
#ifndef FIELDS_H_INCLUDED
#define FIELDS_H_INCLUDED 1

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


class Fields
{
public:
    /**
     * \brief Getting Magnetic field
     * \param v Location of query in cartesian coordinates
     * \return a Vector containing (Bx, By, Bz)
     * \details Must be implemented by children classes
     */
    virtual Vector getB(Vector& v);
    
    /**
     * \brief Whether the location v is at the limiter
     * \note Must be inherited by child class
     */
    virtual bool isLimiter(Vector v);
    
    /**
     * \brief writes the values of matrix components and Psi along row w
     */
    void writeMat(std::string& filename, MatDoub& m);
    
    void writeVec(std::string& filename, VecDoub& v);
    
protected:
    friend class FieldsTests;
};

#endif
