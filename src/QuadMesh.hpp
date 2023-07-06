/**
 * \file    QuadMesh.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    September, 2020
 *
 * \brief   Declares the QuadMesh class
 * \details Inherited from Mesh class
 */

#ifndef quadmesh_h
#define quadmesh_h

#include "Mesh.hpp"
#include "TriagMesh.hpp"

struct Quad{
    MeshPoint ll_, lr_, ul_, ur_;
    LineSeg left, top, right, bottom;
    
    /**
     * \brief Constructor for Quadrangles
     * \note constructs the line segments for the 4 edges.
     */
    Quad(MeshPoint ll, MeshPoint lr, MeshPoint ul, MeshPoint ur);
    
    /**
     *\brief Return the centroid of the quad
     */
    MeshPoint centroid();
    
    /**
     *\brief Whether the quad contains point p
     *\param p the MeshPoint in question
     *\param intersect if point is outside, modified to indicate which side of the
     *      quad the point is next to: {l, t, r, b}
     *\return true if p is in the quad
     */
    bool contains(MeshPoint p, int& intersect);
    
    /**
     * \brief Define how a vector is printed
     */
    friend std::ostream& operator<<(std::ostream &os, const Quad& q);
};

class QuadMesh : public Mesh
{
public:
    /**
     *\brief Constructor for Quadrangle grid
     *\param crx x coordinates of grid points {lower left, lower right, upper left, upper right}
     *\param cry y coordinates of grid points {lower left, lower right, upper left, upper right}
     *\param nbx index of the x coordinate of neighbor cells {left, top, right, bottom}
     *\param nby index of the y coordinate of neighbor cells {left, top, right, bottom}
     *\param vals values defined on the grid
     */
    QuadMesh(std::vector<MatDoub>& crx, std::vector<MatDoub>& cry,
             std::vector<MatDoub>& nbx, std::vector<MatDoub>& nby,
    		std::vector<MatDoub>& vals);

    /**
     *\brief Get number of vals per row
     */
    double nx();
    
    /**
     *\brief Get number of vals per col
     */
    double ny();

    /**
     * \brief Whether the mesh has been initialized
     */
    bool isRaw();
    
    /**
     *\brief Look for which quadrangle the point is in by "walking" through the mesh
     *\param p The point to search for
     *\param init Index of the quad to start with (iw, ih); updated to the search result
     *\return the (ix, iy) index history of the result quadrangle
     */
    std::vector<size_t> search(MeshPoint p, std::vector<size_t>& init);


    /**
     *\brief Step-function interpolation on the quad grid
     *\param p        Query point for interpolation
     *\param init Index of the quad to start with (iw, ih); updated to the search result
     *\return A vector containing each of the values defined on the grid.
     *\details Search for which quad the point is in, and returns the val_ in the quad
     */
    VecDoub interp(MeshPoint p, std::vector<size_t>& init);


    /**
     *\brief Prints the grid and corresponding values
     */
    void print(std::string& fname);
    
    /**
     *\brief Construct a Quad structure for the given q index
     */
    Quad ConstructQuad(std::vector<size_t>& qNow);

private:
    // coordinate of vertices
	std::vector<MatDoub> crx_;
	std::vector<MatDoub> cry_;
    
    // index of physical neighbors
    std::vector<MatDoub> nbx_; // {l, t, r, b}
    std::vector<MatDoub> nby_;
    
    // values defined in each cell
	std::vector<MatDoub> vals_;
    
    friend class MeshTests;
    
    int nx_, ny_;
    
	bool initialized_;
    
    /**
     *\brief Walk to the next quad if p is not in current quad
     *\param p The point to search for
     *\param qNow The (ix, iy) index of the current quadrangle. If p is not in the current quad,
     *          value is modified to the next quad to search
     *\return true if p is in qNow (search is complete)
     */
    bool walk(MeshPoint p, std::vector<size_t>& qNow);
    

};

#endif // quadmesh included
