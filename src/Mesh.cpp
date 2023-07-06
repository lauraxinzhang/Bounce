/**
 * \file    Mesh.cpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    September, 2020
 *
 * \brief   Implements the Mesh class
 * \details 
 */

#include "Mesh.hpp"

std::ostream& operator<<(std::ostream &os, const MeshPoint& p)
{
    os << p.x_ << ", " << p.y_ << std::endl;
    return os;
}

LineSeg::LineSeg(MeshPoint& pa, MeshPoint& pb)
        :pa_(pa), pb_(pb) // using copy constructor
{
    // nothing else to do
}

bool LineSeg::parallel(LineSeg &line_b)
{
//    MeshPoint diff_a, diff_b;
    MeshPoint diff_a = pa_ - pb_;
    MeshPoint diff_b = line_b.pa_ - line_b.pb_;
    double result = diff_a.cross(diff_b);
    return (fabs(result) < 1e-15);
}

double LineSeg::signedArea(MeshPoint& pc)
{
    MeshPoint diff_ba, diff_ca;
    diff_ba = pb_ - pa_;
    diff_ca = pc - pa_;
    return diff_ba.cross(diff_ca) / 2;
}

MeshPoint LineSeg::intersect(LineSeg &line_b)
{
    bool par = parallel(line_b);
    MeshPoint result;
    if (par){
        double nan = std::nan("0");
        MeshPoint rtn(nan, nan);
        result = rtn;
    }
    else{
//        MeshPoint diff_a, diff_b;
        MeshPoint diff_a = pb_ - pa_;
        MeshPoint diff_b = line_b.pa_ - line_b.pb_;
        double denom = 1 / diff_a.cross(diff_b);
        
//        MeshPoint diff_heads;
        MeshPoint diff_heads = pa_ - line_b.pa_;
        double s = diff_b.cross(diff_heads) * denom;
        double t = diff_heads.cross(diff_a) * denom;
        MeshPoint rtn(s, t);
        result = rtn;
    }
    return result;
}

MeshPoint LineSeg::intersectCoord(LineSeg &line_b)
{
    MeshPoint bary = intersect(line_b);
    MeshPoint diff = (pb_ - pa_) * bary.x_;
    MeshPoint rtn = pa_ + diff;
    return rtn;
}

bool LineSeg::isCross(LineSeg &line_b)
{
    MeshPoint inter = intersect(line_b);
//    bool head = inter.x_>= 0 && inter.x_ <= 1;
//    bool tail = inter.y_>= 0 && inter.y_ <= 1;
    bool head = inter.x_>= PRECISION;
    bool head2 = inter.x_- 1 <= PRECISION;
    bool tail = inter.y_>= PRECISION;
    bool tail2 = inter.y_ - 1 <= PRECISION;
    return head && head2 && tail && tail2;
}

bool LineSeg::isCrossEx(LineSeg &line_b)
{
    MeshPoint inter = intersect(line_b);
//    bool head = inter.x_>= 0 && inter.x_ <= 1;
//    bool tail = inter.y_>= 0 && inter.y_ <= 1;
    bool head = inter.x_> PRECISION;
    bool head2 = inter.x_- 1 < PRECISION;
    bool tail = inter.y_> PRECISION;
    bool tail2 = inter.y_ - 1 < PRECISION;
    return head && head2 && tail && tail2;
}

double LineSeg::length()
{
    MeshPoint diff = pa_ - pb_;
    double rsquared = diff.x_ * diff.x_ + diff.y_ * diff.y_;
    return sqrt(rsquared);
}

Mesh::Mesh()
    : initialized_(false) {}

// double Mesh::size()
// {
//     return 0;
// }

//std::vector<size_t> Mesh::search(MeshPoint p, size_t init)
//{
//    try {
//        throw MeshException(6);
//    } catch (MeshException& e) {
//        std::cerr << "in Mesh::search" << std::endl;
//    }
//    std::vector<size_t> rtn{init};
//    return rtn;
//}

double Mesh::interp(MeshPoint p, size_t& init)
{
    try {
        throw MeshException(6);
    } catch (MeshException& e) {
        std::cerr << "in Mesh::interp" << std::endl;
    }
//    VecDoub rtn = {0};
    return 0;
}

VecDoub Mesh::interp(MeshPoint p, std::vector<size_t>& init)
{
    try {
        throw MeshException(6);
    } catch (MeshException& e) {
        std::cerr << "in Mesh::interp" << std::endl;
    }
    VecDoub rtn = {0};
    return rtn;
}

void Mesh::print(std::string& fname)
{
    try {
        throw MeshException(6);
    } catch (MeshException& e) {
        std::cerr << "in Mesh::print" << std::endl;
    }
}

bool Mesh::isRaw()
{
    std::cerr << "Warning: isRaw is being called in Mesh base class." << std::endl;
    return false;
}


