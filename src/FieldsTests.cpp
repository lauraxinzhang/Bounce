/**
 * \file    FieldsTests.cpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2020
 *
 * \brief   Implements the FieldsTests class
 * 
 */
#include "FieldsTests.hpp"

simpleFields::simpleFields(int nw, int nh)
        :GEquilibrium(nw, nh) {}

void simpleFields::setLinear(double rdim, double zdim, double kr, double kz)
{
    int nw = gfile_.nw;
    int nh = gfile_.nh;

    gfile_.rdim = rdim;
    gfile_.zdim = zdim;
    gfile_.rcentr = rdim / 2;
    gfile_.rleft = 0;
    gfile_.zmid = zdim / 2;
    
    gfile_.sibry = 1;
    gfile_.simag = 0;
    
    double dsi = (gfile_.sibry - gfile_.simag)/ (gfile_.nw - 1);

    for (int i = 0; i < nw; ++i){
        gfile_.siGrid_.at(i) = gfile_.simag + i * dsi;
        gfile_.fpol_.at(i) = gfile_.siGrid_.at(i); // * gfile_.siGrid_.at(i);
    }

    double dr = gfile_.rdim / (gfile_.nw - 1);
    // std::cout << "dr " << dr << std::endl;
    for (int i = 0; i < nw; ++i){     // 'column' index
        gfile_.rGrid_.at(i) = gfile_.rleft + i * dr;
    }

    double dz = gfile_.zdim / (gfile_.nh - 1);
    // std::cout << "dz " << dz << std::endl;
    for (int j = 0; j < nh; ++j){ // 'row' index
       gfile_.zGrid_.at(j) = (gfile_.zmid - gfile_.zdim/2.0) + j * dz;
    }
    
    for (int i = 0; i < nh; ++i){    // 'column' index, in r
       for (int j = 0; j < nw; ++j){ // 'row' index, in h
           double rnow = gfile_.rGrid_.at(j);
           double znow = gfile_.zGrid_.at(i);
           gfile_.psirz_.at(j, i) = kr * rnow + kz * znow;
       }
    }
}

void FieldsTests::runAll()
{
    testLinear();
    testGrid();
    testGamma();
    testDMat();
}

void FieldsTests::testLinear()
{
    bool grid = testGrid();
    if ( ! grid) std::cout << "Grid are not setup correctly" << std::endl;
    
    bool psiRZ = testPsiRZ();
    if (! psiRZ) std::cout << "psiRZ is not calculated correctly" << std::endl;
}

bool FieldsTests::testGrid()
{
    simpleFields linear(5, 6);
    linear.setLinear(1, 1, 1, 1);

    VecDoub rgridExp = {0, 0.25, 0.5, 0.75, 1};
    VecDoub zgridExp = {0, 0.2, 0.4, 0.6, 0.8, 1};

    VecDoub sigridExp = rgridExp;

    bool rgEq = vecEqual(rgridExp, linear.gfile_.rGrid_);
    bool zgEq = vecEqual(zgridExp, linear.gfile_.zGrid_);
    bool sigEq = vecEqual(sigridExp, linear.gfile_.siGrid_);
    bool fpolEq = vecEqual(sigridExp, linear.gfile_.fpol_);

    bool rtn = (rgEq && zgEq && sigEq && fpolEq);
    std::string what("LinearField grid");
    result(rtn, what);
    if (!rtn) {
        if (! rgEq) {
            std::cout << "rGrid is wrong" << std::endl;
            std::cout << linear.gfile_.rGrid_ << std::endl;
        }
        if (! zgEq) {
            std::cout << "zGrid is wrong. Got: " << std::endl;
            std::cout << linear.gfile_.rGrid_;
            std::cout << std::endl;
        }
        if (! sigEq) {
            std::cout << "siGrid is wrong" << std::endl;
            std::cout << linear.gfile_.siGrid_;
            std::cout << std::endl;
        }
        if (! fpolEq) {
            std::cout << "fpol is wrong" << std::endl;
            std::cout << linear.gfile_.fpol_;
            std::cout << std::endl;
        }
    }
    return rtn;
}

bool FieldsTests::testPsiRZ()
{
    simpleFields linear(5, 6);
    linear.setLinear(1, 1, 1, 2);
    bool rtn(1);
    for (int iw = 0; iw < linear.gfile_.nw; ++iw){
        for (int ih = 0; ih < linear.gfile_.nh; ++ih){
            double psiNow = linear.gfile_.psirz_.at(iw, ih);
            double rnow = linear.gfile_.rGrid_.at(iw);
            double znow = linear.gfile_.zGrid_.at(ih);
            double psiExp = rnow + 2 * znow;
            rtn = rtn && doubEqual(psiNow, psiExp);
        }
    }
    if (rtn) std::cout << "PsiRZ is setup correctly in linear field" << std::endl;
    return rtn;
}

bool FieldsTests::testGamma()
{
    simpleFields linear(5, 6);
    int zero(0);
    int odd(3);
    int even(4);
    
    bool rtn(1);
    rtn = rtn && linear.gamma(zero) && linear.gamma(even);
    rtn = rtn && (-1 * linear.gamma(odd));
    if (rtn) std::cout << "Fields.gamma is behaving correctly." << std::endl;
    return rtn;
}


bool FieldsTests::testDMat()
{
    simpleFields linear(5, 6);
    
    VecDoub grid = {0, 0.25, 0.5, 0.75, 1};
    MatDoub D = linear.Dij(grid);
    
    VecDoub row1 = { -17./3,  8.,   -4.,  8./3, -1.   };
    VecDoub row2 = { -2.,    -2./3,  4., -2.,    2./3 };
    VecDoub row3 = {  1.,   -4.,    0.,  4.,   -1.   };
    VecDoub row4 = { -2./3,  2.,   -4.,  2./3,  2.   };
    VecDoub row5 = {  1.,   -8./3,  4., -8.,    17./3 };
    
//    VecDoub row1 = { -17./3, -2.,    1.,  -2./3,  1    };
//    VecDoub row2 = {  8.,    -2./3, -4.,  2.,    -8./3 };
//    VecDoub row3 = { -4.,     4.,    0.,  -4.,    4.   };
//    VecDoub row4 = {  8./3,  -2.,    4.,  2./3,  -8.   };
//    VecDoub row5 = { -1.,     2./3, -1.,  2.,    17./3 };
    
    std::vector<VecDoub> matExp = {row1, row2, row3, row4, row5};
    
    bool rtn(1);
    for (int i = 0; i < D.getH(); ++i){
        VecDoub rowNow = D.getRow(i);
        rtn = rtn && vecEqual(matExp.at(i), rowNow);
    }
    if (rtn) {
        std::cout << "D matrix is calculated correctly" << std::endl;
    } else {
        std::cout << "D matrix is wrong. got: " << std::endl;
        std::cout << D << "expected: " << std::endl;
        std::cout << row1 << row2 << row3 << row4 << row5;
        
    }
    return rtn;
}

bool FieldsTests::testB()
{
//    simpleFields linear(5, 6);
//    double kr(1), kz(2);
//    linear.setLinear(1, 1, kr, kz); // val = r + z;
    
    return false;
}
