/**
 * \file    UniMesh.cpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    October, 2020
 *
 * \brief   Implements the UniMesh class; A class inherited from
 *          Mesh, specifically for uniform grids
 * \details 
 */

#include "UniMesh.hpp"

UniMesh::UniMesh(VecDoub& grid, VecDoub& val)
		:initialized_(true), //reversed_(false),
         xgrid_(grid), ygrid_(0), val_(val), dy_(0), dim_(1), linear_(USELINEAR)
{
    dx_ = xgrid_.at(1) - xgrid_.at(0);
}

UniMesh::UniMesh(VecDoub& xgrid, VecDoub& ygrid, MatDoub& val2D)
        :initialized_(true), //reversed_(false),
         xgrid_(xgrid), ygrid_(ygrid),
         val_(val2D), dim_(2), linear_(USELINEAR)//, sprows_(ygrid.size())
{
    dx_ = xgrid_.at(1) - xgrid_.at(0);
    dy_ = ygrid_.at(1) - ygrid_.at(0);
//    for (int i = 0; i < sprows_.size(); i++){
//        // ensure that these pointers are null instead of undefined behavior
//        sprows_.at(i) = nullptr;
//    }
}

int UniMesh::nx()
{
    return val_.getW();
}

int UniMesh::ny()
{
	return val_.getH();
}

bool UniMesh::isRaw()
{
	return !initialized_;
}

double UniMesh::interp(MeshPoint p, size_t init)
{
    if (linear_){
//        std::cerr << "linear interp" << std::endl;
        if (dim_ == 1){
            return interp1DLin( p.x_, xgrid_, val_.getRow(0));
        } else {
            return interp2DLin(p);
        }
    } else {
        std::cerr << "BS interp disabled." << std::endl;
        exit(0);
//        if (dim_ == 1){
//            return interp1DBS( p.x_, xgrid_, val_.getRow(0));
//        } else {
//            return interp2DBS(p);
//        }
    }
}

int UniMesh::search(double query, VecDoub &grid)
{
    double h = grid.at(1) - grid.at(0); // negative if descending
    double numspaces = ( query - grid.front() ) / h; // always positive
    return floor(numspaces);
}

double UniMesh::interp2DLin(MeshPoint p)
{
//    int iy = floor((p.y_ - ygrid_.front()) / dy_) ;
    int iy = search(p.y_, ygrid_);
    
//	std::cerr << p.y_ << ',' << ygrid_.back() << std::endl;
    double rtn(0);
    try {
        double y0 = ygrid_.at(iy);
//        std::cerr << ygrid_.size() << ',' << iy << std::endl;
        double y1 = ygrid_.at(iy + 1);
        double y_val_0 = interp1DLin(p.x_, xgrid_, val_.getRow(iy));
        double y_val_1 = interp1DLin(p.x_, xgrid_, val_.getRow(iy + 1));
        
        VecDoub ygrid = {y0, y1};
        VecDoub valgrid = {y_val_0, y_val_1};
        
//        std::cerr << "in interp2D. query: " << p.y_ << std::endl;
//    	std::cerr << "grid: " << ygrid << std::endl; 
        rtn = interp1DLin(p.y_, ygrid, valgrid);
        
   } catch (std::out_of_range& oor) {
	std::cerr << "std oor caught" << std::endl;
        if (fabs(p.y_ - ygrid_.back()) < PRECISION){
            rtn = interp1DLin(p.x_, xgrid_, val_.getRow(iy));
        } else {
//            std::cerr << oor.what() << std::endl;
            throw MeshException(9);
//            e.what();
        }
    }
    return rtn;
}

double UniMesh::interp1DLin(double query, VecDoub& grid, VecDoub& val)
{
    
//    int i = floor( ( query - grid.front() ) / h) ;
    int i = search(query, grid);
    double result(0);
    //std::cerr << "in interp1D. found i: " << i << std::endl;
    try {
        double x0 = grid.at(i);
        double f0 = val.at(i);
        double x1 = grid.at(i+1);
        double f1 = val.at(i+1);
        
        double numerator = f0 * (x1 - query) + f1 * ( query - x0 );
        double h = x1 - x0;
        result = numerator / h;
        
    } catch (std::out_of_range& oor) {
        if (fabs(query - grid.back()) < PRECISION){
            result = val.back();
//            return result;
        } else {
//            std::cerr << oor.what() << std::endl;
            std::cerr << "in interp1D. query: " << query << std::endl;
            std::cerr << "grid: " << grid << std::endl; 
            throw MeshException(9);


//            e.what();
        }
    }
    return result;
}


void UniMesh::print(std::string& fname)
{
	return;
}

int XYMesh::search(double query, VecDoub& grid){
    int i;
    if (grid.at(1) > grid.at(0)){
        auto lower = std::lower_bound(grid.begin(), grid.end(), query);
        i = std::distance(grid.begin(), lower) - 1;
    } else {
        // the grid is in descending order. use reverse iterators instead
        auto lower = std::lower_bound(grid.rbegin(), grid.rend(), query);
        i = std::distance( lower, grid.rend()) - 1;
    }
    return i;
}
