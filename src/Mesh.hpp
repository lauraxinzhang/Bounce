/**
 * \file    Mesh.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    September, 2020
 *
 * \brief   Declares the MeshPoint & LineSeg structures & 
 *          the Mesh virtual interface
 * \details 
 */

#ifndef mesh_h
#define mesh_h

#include "delaunator.hpp"
#include "Containers.hpp"
#include <cmath>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <stdio.h>
#include <algorithm>

typedef std::vector<double> VecDoub;
struct MeshPoint
{
    double x_;
    double y_;
    
    /**
     * \brief Default constructor
     */
    inline MeshPoint():x_(0), y_(0) {}
    
    /**
     * \brief Constructor
     * \param x x coordinate
     * \param y y coordinate
     */
    inline MeshPoint(const double& x, const double& y):x_(x), y_(y){}
    
    /**
     * \brief Copy Constructor
     */
    inline MeshPoint(const MeshPoint& p): x_(p.x_), y_(p.y_) {}
    
    /**
     * \brief Comparison operator
     */
    inline bool operator==(MeshPoint& p)
    {
        return p.x_ == x_ && p.y_ == y_;
    }
    
    /**
     * \brief Subtraction 
     * \details 2D vector subtraction
     * \return (x2 - x1, y2 - y1)
     */
    inline MeshPoint operator-(MeshPoint& pb)
    {
        double x = x_ - pb.x_;
        double y = y_ - pb.y_;
        MeshPoint p(x, y);
        return p;
    }
    
    /**
     * \brief Addition
     * \details 2D vector addition
     * \return (x2 + x1, y2 + y1)
     */
    inline MeshPoint operator+(MeshPoint& pb)
    {
        double x = x_ + pb.x_;
        double y = y_ + pb.y_;
        MeshPoint p(x, y);
        return p;
    }
    
    /**
     * \brief Scalar multiplication
     * \details 2D vector addition
     * \return (kx, ky)
     */
    inline MeshPoint operator*(const double& k)
    {
        double x = k * x_;
        double y = k * y_;
        MeshPoint p(x, y);
        return p;
    }
    
    /**
     * \brief Two-dimensional cross product
     */
    inline double cross(MeshPoint& pb)
    {
        double rtn = x_ * pb.y_ - pb.x_ * y_;
        return rtn;
    }
    
    /**
     * \brief Define how a MeshPoint is printed
     */
    friend std::ostream& operator<<(std::ostream &os, const MeshPoint& p);
};



struct LineSeg
{
    MeshPoint pa_;
    MeshPoint pb_;
    
    /**
     * \brief Constructor
     */
    LineSeg(MeshPoint& pa, MeshPoint& pb);
    
    /**
     * \brief Check whether two line segments are parallel
     */
    bool parallel(LineSeg& line_b);
    
    /**
     *\brief Value of signed area of the triangle formed by the current linesegment and point c
     *\param pc the 3rd MeshPoint to form triangle with
     *\note Return value is positive if pc is "left-of" line segment
     */
    double signedArea(MeshPoint& pc);
    
    /**
     *\brief Computes the intersection of two lines.
     *\param line_b Query line segment.
     *\return MeshPoint(s, t) of the intersection, in the basis of the line segments. See nr3 ch21.4.
     */
    MeshPoint intersect(LineSeg& line_b);
    
    /**
     *\brief Intersect of two line segments in their own coordinates
     */
    MeshPoint intersectCoord(LineSeg& line_b);
    
    /**
     *\brief Whether the two line SEGMENTS intersect
     *\param line_b Query line segment
     *\return true if the computed intersect lies *on* the line segment, including end points.
     */
    bool isCross(LineSeg& line_b);
    
    /**
     *\brief Whether the two line SEGMENTS intersect
     *\param line_b Query line segment
     *\return true if the computed intersect lies *on* the line segment, EXCLUDING end points.
     */
    bool isCrossEx(LineSeg& line_b);
    
    /**
     *\brief Length of the line segment
     */
    double length();
};



/**
 *\brief Housekeeping. Allows for exiting without memory leaks
 */
struct MeshException : public std::exception {
    int code_; ///< Exit code
    MeshException(int code): code_(code) {} // trivial constructor
    
    const char * what () const throw () {
        std::cerr << "MeshException encountered: " << std::endl; 
        switch (code_){
            case 0: std::cerr << "Fatal: Mesh not initialized" << std::endl;
                exit(0);
            case 1:
                std::cerr << "Fatal: Illegal parameters provided in TriagMesh contructor. Coordinates and function values dimension mismatch." << std::endl;
                std::cerr << "Exiting..." << std::endl;
                exit(1);
            case 2:
                std::cerr << "Fatal: Point is not in triangle but no intersection found with current triangle." << std::endl;
                std::cerr << "Exiting..." << std::endl;
                exit(2);
            case 3:
//                std::cerr << "Fatal: Search target point is outside domain. Exiting...";
                std::cerr << "Search target is outside trianglation." << std::endl;
                std::cerr << std::endl;
                exit(3);
            case 4:
                std::cerr << "Fatal: Triangle index out of range. Exiting..." << std::endl;
                std::cerr << "Called by Mesh::edgesOfTriag" << std::endl;
                exit(4);
            case 5:
                std::cerr << "Fatal: Edge index out of range. Exiting..." << std::endl;
                std::cerr << "Called by Mesh::triagOfEdge" << std::endl;
                exit(5);
            case 6:
                std::cerr << "Warning: A virtual function in the Mesh base class is called." << std::endl;
//                exit(6);
            case 7:
                std::cerr << "In walk: Domain boundry encoutered on QuadMesh" << std::endl;
            case 8:
                std::cerr << "Fatal: in ConstructQuad. requested quad does not exist." << std::endl;
                exit(8);
            case 9:
                std::cerr << "Error encountered in Uniform Mesh. ";
                std::cerr<< "The query point for interpolation is outside of provided domain" << std::endl;
                exit(9);
            case 10:
                std::cerr << "Veroni cells for boundary points are not implemented." << std::endl;
                exit(10);
        }
        return "Uncaught exceptions";
   }
};

class Mesh
{
public:

    /**
     * \brief Default constructor
     * \details Initializes initialized_ to false
     */
    Mesh();

    /**
     *\brief Get number of coordinate pairs
     */
    // double size();


//    /**
//     *\brief Look for which triangle the point is in by "walking" through the triangulation
//     *\param p The point to search for
//     *\param init Index of the triangle to start with
//     */
//    virtual std::vector<size_t> search(MeshPoint p, size_t init);

    /**
     *\brief Linear interpolation on a grid
     *\param p        Query point for interpolation
     *\param init Index of the initial angle to start searching in
     *\details Search for which cell the point is in, then interpolate
     *TODO: this signature is inconsistent with quad. change to allow for vectors
     *here and in Triag.
     */
    virtual VecDoub interp(MeshPoint p, std::vector<size_t>& init);
    
    /**
     *\brief Linear interpolation on triangular grid
     *\param p        Query point for interpolation
     *\param init Index of the initial triangle to start searching in
     *\details Search for which cell the point is in, then interpolate
     *TODO: this signature is inconsistent with quad. change to allow for vectors
     *here and in Triag.
     */
    virtual double interp(MeshPoint p, size_t& init);

    virtual void print(std::string& fname);

    /**
     * \brief Whether the mesh has been initialized
     */
    virtual bool isRaw();

protected:
    // TODO: do these need to be in the parent class?
    bool initialized_;
    
    friend class MeshTests;

};


#endif /* mesh_h */
