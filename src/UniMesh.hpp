/**
 * \file    UniMesh.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    October, 2020
 *
 * \brief   Declares the UniMesh class; A class inherited from
 *          Mesh, specifically for uniform grids
 * \details 
 */

#ifndef unimesh_h
#define unimesh_h

#include "Constants.hpp"
#include "Mesh.hpp"
#include "TriagMesh.hpp"
//#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

//typedef boost::math::interpolators::cardinal_cubic_b_spline<double> bSpline;

class UniMesh : public Mesh
{
public: 

	UniMesh(VecDoub& grid, VecDoub& val);
    
    UniMesh(VecDoub& xgrid, VecDoub& ygrid, MatDoub& val2D);

	int nx();

	int ny();

	bool isRaw();

	double interp(MeshPoint p, size_t init = 0);
    
    double interp2DLin(MeshPoint p);

    double interp1DLin(double query, VecDoub& grid, VecDoub& val);
    
//    /**
//     * \brief 2D interpolation with boost library BSpline
//     */
//    double interp2DBS(MeshPoint p);
//    
//    void prep2DSpline();
//    
//    /**
//     * \brief 1D interpolation with boost library BSpline
//     */
//    double interp1DBS(double query, VecDoub& grid, VecDoub& val);
    
    virtual int search(double query, VecDoub& grid);

	void print(std::string& fname);

protected:
	bool initialized_;

    
    VecDoub xgrid_;
    VecDoub ygrid_;
    MatDoub val_;
    
    int dim_; // dimension of the grid. 1D or 2D supported.
    
    double dx_;
    double dy_;
    bool linear_;
//    std::vector<bSpline*> sprows_;

    friend class MeshTests;
};

class XYMesh : public UniMesh
{
public:
    inline XYMesh(VecDoub& grid, VecDoub& val)
                : UniMesh(grid, val)
    {
        linear_ = true; // only linear interpolation. BSpline does not work
    }
    
    inline XYMesh(VecDoub& xgrid, VecDoub& ygrid, MatDoub& val2D)
                : UniMesh(xgrid, ygrid, val2D)
    {
        linear_ = true;
    }
    
    int search(double query, VecDoub& grid);
protected:
    bool reversed_; // whether the support grids are ins descending order

};

#endif
