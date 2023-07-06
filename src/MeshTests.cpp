/**
 * \file    MeshTests.cpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    October, 2020
 *
 * \brief   Implements the MeshTests class
 */
#include "MeshTests.hpp"

MeshTests::MeshTests()
        : uniMesh_(nullptr), xgrid_(nullptr), ygrid_(nullptr),
        val1D_(nullptr), val2D_(nullptr), nx_(0), ny_(0)
{}

void MeshTests::runALL()
{
    testUniMesh1D();
    testUniMesh2D();
    testUniInterp();
    testTriagInterp();
}

void MeshTests::prepUniMesh1D()
{
    if (uniMesh_){
        delete uniMesh_;
    }
    nx_ = 6;
    ny_ = 1;
    VecDoub xx = linspace(-3.7, 0, nx_);
    VecDoub val = linspace(0, 10, nx_);
    
    for (int ix = 0; ix < nx_; ix++){
        val.at(ix)  = 2 * xx.at(ix);
    }
    xgrid_ = new VecDoub (xx);
    val1D_ = new VecDoub (val);
    uniMesh_ = new UniMesh((*xgrid_), (*val1D_));
    
    return;
}

void MeshTests::prepUniMesh2D()
{
    if (uniMesh_){
        delete uniMesh_;
    }
    nx_ = 6;
    ny_ = 15;
    
    VecDoub xx = linspace(0, 5, nx_);
    VecDoub yy = linspace(0, 9, ny_);
    
    
    MatDoub val(nx_, ny_);
    for (int ix = 0; ix < nx_; ++ix){
        for (int iy = 0; iy < ny_; iy++){
            val.at(ix, iy)  = xx.at(ix) - 2 * yy.at(iy);
        }
    }
    
    uniMesh_ = new UniMesh(xx, yy, val);
    xgrid_ = new VecDoub(xx);
    ygrid_ = new VecDoub(yy);
    val2D_ = new MatDoub(val);
    return;
}

bool MeshTests::testUniMesh1D()
{
    if (!uniMesh_ || uniMesh_-> ny() != 1){
        prepUniMesh1D();
    }
    std::string what("UniMesh1D - Constructor");
    
    bool size = (uniMesh_ -> nx() == nx_) && (uniMesh_ -> ny() == ny_);
    
    bool gridEq = vecEqual( (*xgrid_), (*uniMesh_).xgrid_ );
    gridEq = gridEq && vecEqual( (*val1D_), (*uniMesh_).val_.getRow(0) );
    bool ans = size && gridEq;
    result(ans, what);
    return ans;
}

bool MeshTests::testUniMesh2D()
{
    if (!uniMesh_ || uniMesh_-> ny() == 1 || !ygrid_){
        prepUniMesh2D();
    }
    std::string what("UniMesh2D - Constructor");
    
    
    bool size = (uniMesh_ -> nx() == nx_) && (uniMesh_ -> ny() == ny_);
    bool gridEq = vecEqual( (*xgrid_), (*uniMesh_).xgrid_ );
    gridEq = gridEq && vecEqual( (*ygrid_), (*uniMesh_).ygrid_ );
    bool valEq = matEqual( (*val2D_), (*uniMesh_).val_ );
    bool ans = size && gridEq && valEq;
    result(ans, what);
    return ans;
}

bool MeshTests::testUniInterp()
{
    std::string what1D("UniMesh - 1D Interpolation");
    std::string what2D("UniMesh - 2D Interpolation");
    
    // first test 1D intepolation
    if (!uniMesh_ || uniMesh_-> ny() != 1){
        prepUniMesh1D();
    }
    bool interp1D(true);
    
    //test interpolation on random points
    VecDoub q_x = randArray(xgrid_->front(), xgrid_->back(), 10);
    for (int i = 0; i < q_x.size(); ++i){
        MeshPoint q(q_x.at(i), 0);
        double out = uniMesh_ -> interp(q);
        interp1D = interp1D && doubEqual(out, 2 * q_x.at(i));
    }
    // test interpolation on grid points
    for (int i = 0; i < xgrid_->size(); ++i){
        double x = xgrid_ -> at(i);
        MeshPoint q( x, 0);
        double out = uniMesh_ -> interp(q);
        interp1D = interp1D && doubEqual(out, 2 * x);
    }
    
    result(interp1D, what1D);
    
    prepUniMesh2D();
    bool interp2D(true);
    q_x = randArray(xgrid_->front(), xgrid_->back(), 10);
    VecDoub q_y = randArray(ygrid_->front(), ygrid_->back(), 10);
    for (int ix = 0; ix < q_x.size(); ++ix){
        for (int iy = 0; iy < q_y.size(); iy++){
            double x = q_x.at(ix);
            double y = q_y.at(iy);
            MeshPoint p(x, y);
            double exp = x - 2 * y;
            double out = uniMesh_ -> interp(p);
            
            interp2D = interp2D && doubEqual(exp, out);
        }
    }
    for (int ix = 0; ix < xgrid_->size(); ++ix){
        for (int iy = 0; iy < ygrid_->size(); iy++){
            double x = xgrid_->at(ix);
            double y = ygrid_->at(iy);
            MeshPoint p(x, y);
            double exp = x - 2 * y;
            double out = uniMesh_ -> interp(p);
            
            interp2D = interp2D && doubEqual(exp, out);
        }
    }
    result(interp2D, what2D);
    
    return interp1D && interp2D;
}

bool MeshTests::testTriagInterp()
{
    std::string whatTriag("TriagMesh - Interpolation");
    VecDoub coords = {-1, 1,  1, 1,   1, -1,   -1, -1,
                    0, 0,
                    -1, 0,
                    0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -0.5, 0.5};
    std::vector<int> mask = {2, 2, 2, 2, 0, 2, 0, 0, 0, 0};
    // build value array
    VecDoub val(coords.size()/2);
    double kx(1), ky(1.6);
    for (int ix = 0; ix < coords.size(); ix += 2){
        double x = coords.at(ix);
        double y = coords.at(ix + 1);
        double z = kx * x + ky * y;
        val.at(ix / 2) = z;
    }
    TriagMesh mesh(coords, val);
    
    VecDoub random = randArray(-1, 1, 20);
    
    VecDoub q_x(random.begin(), random.begin()+10);
    VecDoub q_y(random.begin() + 10, random.end());
    
    bool interpTriag(true);
    size_t init(0);
    // test random points
    for (int ix = 0; ix < q_x.size(); ++ix){
        for (int iy = 0; iy < q_y.size(); iy++){
            double x = q_x.at(ix);
            double y = q_y.at(iy);
            MeshPoint p(x, y);
            double exp = kx * x + ky * y;
            double out = mesh.interp(p, init);
            interpTriag = interpTriag && doubEqual(exp, out);
        }
    }
    
    // test interpolation on grid points
    for (int i = 0; i < coords.size(); i+=2){
        double x = coords.at(i);
        double y = coords.at(i + 1);
        MeshPoint q(x, y);
        size_t init(0);
        double out = mesh.interp(q, init);
        double exp = val.at(i / 2);
        interpTriag = interpTriag && doubEqual(out, exp);
    }
    result(interpTriag, whatTriag);
    
    return interpTriag;
}
