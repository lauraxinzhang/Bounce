/**
 * \file    TriagMesh.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    September, 2020
 *
 * \brief   Declares the TriagMesh class
 * \details Inherited from Mesh class
 */

#ifndef triagmesh_h
#define triagmesh_h

#include "Mesh.hpp"
#include <set>

class TriagMesh : public Mesh
{
public:
    
    /**
     *\brief Constructor for Mesh class. Constructs the triangulation from input coordinates
     *\param coords Set of coordinates for input points, as one vector {x1, y1, x2, y2, ...}
     *\param val Function value to be interpolated.
     *\note The size of the val vector is required to be half of that of coords; only a single set of values is supported at this time.
     */
    TriagMesh(VecDoub& coords, VecDoub& val);
    
    /**
     *\brief Constructor for Mesh class. Constructs the triangulation from input coordinates
     *\param coords Set of coordinates for input points, as one vector {x1, y1, x2, y2, ...}
     *\param val Function value to be interpolated.
     *\note The size of the val vector is required to be half of that of coords; only a single set of values is supported at this time.
     */
    TriagMesh(VecDoub& coords, VecDoub& val, std::vector<int>& contour);
    
    /**
     *\brief Get number of coordinate pairs
     */
    size_t size();
    
    /**
     *\brief Get the number of total triangles
     */
    size_t numTriag();
    
    /**
     *\brief Look for which triangle the point is in by "walking" through the triangulation
     *\param p The point to search for
     *\param init Index of the triangle to start with, updated to the search result
     */
    std::vector<size_t> search(MeshPoint p, size_t& init);
    
    /**
     *\brief Walks to the next closest triangle
     *\param p MeshPoint to locate
     *\param t_now Index of current triangle
     *\return index of adjacent triangle to walk to
     */
    size_t walk(MeshPoint p, size_t t_now);
    
    /**
     *\brief Linear interpolation on triangular grid
     *\param p        Query point for interpolation
     *\param init Index of the initial triangle to start searching in, updated to the search result
     *\details Search for which triangle the point is in, then find the barycentric coordinate of the point. The interpolated value is then sum(lambda_i * val_i).
     */
    double interp(MeshPoint p, size_t& init);
    
    /**
     *\brief An overloaded linear interpolation for values defined on the same mesh
     *\param p query point for interpolation
     *\param val values to be interpolated, defined on the same mesh
     *\param init index of initial guess triangle
     */
    double interp(MeshPoint p, VecDoub& val, size_t& init);
    
    /**
     * \brief print the coordiantes of the triangles to file
     * \param fname name of file to output to
     */
    void print(const char* fname);
    
    /**
     * \brief print the coordiantes of the triangles to file
     * \param fname name of file to output to
     */
    void print(std::string& fname);
    
    /**
     * \brief print the coordiantes of the triangles to file
     * \param fname name of file to output to
     * \param datamask mask of which traingles to print
     */
    void print(std::string& fname, std::vector<int> datamask);
    
    /**
     * \brief Print a different set of data associated with the same triangulation
     * \param data Data to print
     * \param datamask Data stencil (from Grid class)
     */
    void print(std::string& fname, VecDoub& data, std::vector<int> datamask);
    
    /**
     * \brief print the coordiantes of the triangles and the corresponding values to file
     * \param fname name of file to output to
     */
    void printVal(std::string& fname);
    
    void printVal(const char* fname, VecDoub& data);
    
    void printVal(std::string& fname, VecDoub& data);
    
    bool isRaw();
    
    /**
     *\brief Find the indices to the points of triangles in the coord list
     *\param t index of triangle
     *\return indices of points
     */
    std::vector<size_t> pointsOfTriag(size_t t);
    
    void printVeroni(std::string& prefix, std::vector<int> datamask);
    void printTriag(std::string& prefix, std::vector<int>& datamask);
    
    void printTriag(std::string& prefix, std::vector<int>& datamask, std::vector<VecDoub>& values);
    
    void printTriag(std::string& prefix, std::vector<int>& datamask, std::vector<bool>& valid, std::vector<VecDoub>& values);
    /**
     * \brief Find index of all edges going into the starting point
     * \param start Starting t index
     * \param eout Vector to write output to. valid incoming edges.
     * \param ebdry Vector to write output to. boundary edges that have no opposite half edge.
     */
    void edgesAroundPoint(size_t start, std::vector<size_t>& eout);
    
    /**
     *\brief Find veroni cell around a vertex
     *\param start Starting t index
     *\param out Vector to write output to. An list of triangles whose circumsenter forms the veroni cell.
     *\note use std::set for faster search
     */
    void findVeroni(size_t start, std::vector<size_t>& out);
    
private:
    VecDoub coords_;
    VecDoub val_; //length is half of the length of coords; Value on each grid point
    delaunator::Delaunator d_;
    
    friend class MeshTests;
    friend class Grid;
    friend class SOL2P;
    
    /**
     *\brief A flag for whether the mesh is properly initialized. false if only default constructor is called.
     */
    bool initialized_;
    
    /**
     *\brief  Find the half edges coresponding to the input triangle
     *\param  t index of triangle
     *\return Indices of the edges of the t-th triangle
     */
    std::vector<size_t> edgesOfTriag(size_t t);
    
    
    
    /**
     *\brief Find the triangle that the given half edge is in
     *\param e index of the half edge
     *\return Index of the resulting triangle (index t that the triangle starts in)
     */
    size_t triagOfEdge(size_t e);
    
    /**
     *\brief Move to the next half edge in the same triangle
     */
    size_t nextHalfedge(size_t e);
    
    /**
     * \brief Move to the previous half edge in the same triangle
     */
    size_t prevHalfedge(size_t e);
    
    /**
     *\brief Give the coordinates of the vertices of triangle
     *\param t index of triangle
     *\return coordinates of vertices
     */
    std::vector<MeshPoint> coordsOfTriag(size_t t);
    
    /**
     *\brief Give the coordinates of the end points of an half edge (head, tail)
     *\param e index of halfedge
     *\return coordinates of endpoints
     */
    LineSeg edgeToLineSeg(size_t e);
    
    /**
     *\brief Find the opposite halfedge of the given edge, and then the triangle it belongs to
     *\param e index of the half edge in the current triangle.
     *\return index of the neighbor triangle; -1 if there's no neighbor
     */
    size_t neighborTriag(size_t e);
    
    
    /**
     *\brief Calculates the barycentric coordinate of point in triangle
     *\param point coordinate {x, y} of query point
     *\param t          index of triangle
     *\return Barycentric coordinate {lambda1, lambda2, lambda3}; expression from wikipedia
     */
    VecDoub barycentric(MeshPoint point, size_t t);
    
    /**
     *\brief Whether a point is in the current triangle
     *\param point coordinate {x, y} of the query point
     *\param t          index of the current triangle
     *\details Calculates the barycentric coordinate. If all coordinates are between [0, 1],
     *           the point is in the triangle.
     */
    bool isInTriag(MeshPoint point, size_t t);
    
    /**
     *\brief Finds the centroid of the triangle
     *\param t Index of the triangle
     *\note  This is a point that is gauranteed to be *inside* the triangle. Its barycentric coordinate
     *          should always be {1/3, 1/3, 1/3}
     */
    MeshPoint centroid(size_t t);
    
    
    /**
     *\brief Finds the circumcenter of the triangle
     *\param t Index of the triangle
     *\note  This is a point that is gauranteed to be *inside* the triangle. Its barycentric coordinate
     *          should always be {1/3, 1/3, 1/3}
     */
    MeshPoint circumcenter(size_t t);
    

    

};

#endif //triagmesh.h included
