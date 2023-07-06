/**
 * \file    FieldsTests.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2020
 *
 * \brief   Declares the FieldsTests class, a set of tests using SimpleFields
 */
 

#ifndef FIELDSTESTS_H_INCLUDED
#define FIELDSTESTS_H_INCLUDED 1

#include "Fields.hpp"
#include "GEquilibrium.hpp"
#include "Constants.hpp"
#include "TestsUtil.hpp"

class simpleFields: public GEquilibrium
{
public:
    /**
     * \brief Constructor for a test configuration
     * \note  All other inherited constructors are not accessible from outside
     * \param nw Number of horizontal grid points
     * \param nh Number of vertical grid points
     */
    simpleFields(int nw, int nh);

    /**
     * \brief fills the gfile_ with a linear profile of psiRZ
     * \note Only rGrid, zGrid, siGrid, and the corresponding fpol and psirz are filled in.
     * \param rdim dimension of the grid along r
     * \param zdim dimension of the grid along z
     * \param kr slope ws.t. r: val = kr * r + kz * z
     * \param kz slope ws.t. z: val = kr * r + kz * z
     */
    void setLinear(double rdim, double zdim, double kr, double kz);

private:
    friend class FieldsTests;
};


class FieldsTests : public TestsUtil
{
public:
    void runAll();
    
    /**
     * \brief A series of tests to see whether the profiles and fields are initiated as expected
     */
    void testLinear();

    /**
     * \brief Tests whether the grids are setup correctly
     */
    bool testGrid();
    
    /**
     *\brief Tests whether psiRZ Matrix is setup correctly
     */
    bool testPsiRZ();
    
    bool testGamma();
    
    /**
     * \brief Tests whether the differentiation matrix D is calculated correctly
     */
    bool testDMat();

    bool testB();

};

#endif
