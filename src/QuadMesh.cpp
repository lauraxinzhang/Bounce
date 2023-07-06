/**
 * \file    QuadMesh.cpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    September, 2020
 *
 * \brief   Implements the QuadMesh class
 * \details 
 */

#include "QuadMesh.hpp"

Quad::Quad(MeshPoint ll, MeshPoint lr, MeshPoint ul, MeshPoint ur)
    : ll_(ll), lr_(lr), ul_(ul), ur_(ur),
        left(ll, ul), top(ul, ur), right(lr, ur), bottom(ll, lr)
{}

MeshPoint Quad::centroid()
{
    MeshPoint left   = (ll_ + ul_) * 0.5; //center of edges
    MeshPoint right  = (lr_ + ur_) * 0.5;
    MeshPoint top    = (ul_ + ur_) * 0.5;
    MeshPoint bottom = (ll_ + lr_) * 0.5;
            
    LineSeg horizontal(left, right); // medians
    LineSeg vertical(top, bottom);
            
    MeshPoint centroid = horizontal.intersectCoord(vertical);
    
    return centroid;
}

bool Quad::contains(MeshPoint p, int &intersect)
{
    MeshPoint center = centroid();
    LineSeg query(center, p);
    bool rtn(true);
    
    std::vector<LineSeg> edges = {left, top, right, bottom};
    
    for (int i = 0; i < edges.size(); ++i){
        LineSeg e = edges.at(i);
        MeshPoint inter = query.intersect(e);
        if (inter.x_ > 0 && inter.x_ < 1) {
            // intersecting with current edge
            intersect = i;
            rtn = false; // point is outside quad if an intersect is found;
            break;
        }
    }
    return rtn;
}

std::ostream& operator<<(std::ostream &os, const Quad& q)
{
    os << q.ll_ << q.lr_ << q.ul_ << q.ur_;
    return os;
}

QuadMesh::QuadMesh(std::vector<MatDoub>& crx, std::vector<MatDoub>& cry,
                 std::vector<MatDoub>& nbx, std::vector<MatDoub>& nby,
    		     std::vector<MatDoub>& vals)
		: crx_(crx), cry_(cry), nbx_(nbx), nby_(nby),
          vals_(vals), initialized_(true)
{
    nx_ = vals_.at(0).getW();
    ny_ = vals_.at(0).getH();
}

double QuadMesh::nx()
{
    return nx_;
}

double QuadMesh::ny()
{
    return ny_;
}

bool QuadMesh::isRaw()
{
	return !initialized_;
}


void QuadMesh::print(std::string& fname)
{
    return;
}

Quad QuadMesh::ConstructQuad(std::vector<size_t>& qNow)
{
    int ix = qNow.at(0);
    int iy = qNow.at(1);
    if (ix == -1 || ix == nx() || iy == -1 || iy == ny()){
        throw MeshException(8);
    }
    // lower left corner
    MeshPoint ll(crx_.at(0).at(ix, iy), cry_.at(0).at(ix, iy));
    // lower right corner
    MeshPoint lr(crx_.at(1).at(ix, iy), cry_.at(1).at(ix, iy));
    // upper left corner
    MeshPoint ul(crx_.at(2).at(ix, iy), cry_.at(2).at(ix, iy));
    // upper right corner
    MeshPoint ur(crx_.at(3).at(ix, iy), cry_.at(3).at(ix, iy));
    
    Quad quad(ll, lr, ul, ur);
    return quad;
}

VecDoub QuadMesh::interp(MeshPoint p, std::vector<size_t>& init)
{
    std::vector<size_t> hist = search(p, init);
    int ix = init.at(0);
    int iy = init.at(1);
    VecDoub values(vals_.size());
    for (int i = 0; i < values.size(); i++){
        values.at(i) = vals_.at(i).at(ix, iy);
    }
    return values;
}

std::vector<size_t> QuadMesh::search(MeshPoint p, std::vector<size_t>& init)
{
    if (isRaw()){
        throw MeshException(0);
    }
    std::vector<size_t> rtn;
    std::vector<size_t> q_now(init);
    
    rtn.push_back(q_now.at(0));
    rtn.push_back(q_now.at(1));

    while (!walk(p, q_now)){
//        std::cout << q_now.at(0) << ',' << q_now.at(1) << std::endl;
        rtn.push_back(q_now.at(0));
        rtn.push_back(q_now.at(1));
    }
    init = q_now;
    return rtn;
}

bool QuadMesh::walk(MeshPoint p, std::vector<size_t>& qNow)
{
    int ix = qNow.at(0);
    int iy = qNow.at(1);
    if (ix == 0 || ix == nx() || iy == 0 || iy == ny()){
        throw MeshException(7);
        return false;
    }
    Quad quad = ConstructQuad(qNow);
    int inter(-1);
    bool inside = quad.contains(p, inter);
    int ix_next(ix), iy_next(iy);
    if (!inside){
        switch (inter) {
            case 0:
                // move to left neighbor
                ix_next = nbx_.at(0).at(ix, iy);
                iy_next = nby_.at(0).at(ix, iy);
                break;
            case 1:
                // move to top neighbor
                ix_next = nbx_.at(1).at(ix, iy);
                iy_next = nby_.at(1).at(ix, iy);
                break;
            case 2:
                // move to right neighbor
                ix_next = nbx_.at(2).at(ix, iy);
                iy_next = nby_.at(2).at(ix, iy);
                break;
            case 3:
                // move to bottom neighbor
                ix_next = nbx_.at(3).at(ix, iy);
                iy_next = nby_.at(3).at(ix, iy);
                break;
        }
        qNow.at(0) = ix_next;
        qNow.at(1) = iy_next;
    }
    return inside;
}

