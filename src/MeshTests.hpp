/**
 * \file    MeshTests.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2020
 *
 * \brief   Declares the MeshTests class, 
 			a set of tests using for meshes
 */
 

#ifndef MESHTESTS_H_INCLUDED
#define MESHTESTS_H_INCLUDED 1

#include "Fields.hpp"
#include "Constants.hpp"
#include "TestsUtil.hpp"
#include "Mesh.hpp"
#include "TriagMesh.hpp"
#include "QuadMesh.hpp"
#include "UniMesh.hpp"

class MeshTests : public TestsUtil
{
public:
	MeshTests();
    
    void runALL();

    void prepUniMesh1D();
    
    void prepUniMesh2D();
    
    bool testUniMesh1D();
    
    bool testUniMesh2D();
    
    bool testUniInterp();
    
    bool testTriagInterp();
    
private:
    UniMesh * uniMesh_;
    
    VecDoub * xgrid_;
    VecDoub * ygrid_;
    
    VecDoub * val1D_;
    MatDoub * val2D_;
    
    int nx_, ny_;
    
    
};



#endif
