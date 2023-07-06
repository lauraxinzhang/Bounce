/**
 * \file    TriagMesh.cpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    September, 2020
 *
 * \brief   Implements the TriagMesh class
 * \details 
 */

#include "TriagMesh.hpp"

TriagMesh::TriagMesh(VecDoub& coords, VecDoub& val)
    :d_(coords),
     coords_(coords), val_(val), initialized_(true)
{
    // Delaunator is constructed in colon initialization
    if (val_.size() != coords_.size() /2){
        try {
            throw MeshException(1);
        } catch (MeshException& e) {
            std::cerr << "Found " << val.size() << " values for ";
            std::cerr << coords_.size() /2 << " pairs of coordinates." << std::endl;
            std::cout << e.what() << std::endl;
        }
    }
}

TriagMesh::TriagMesh(VecDoub& coords, VecDoub& val, std::vector<int>& contour)
    :d_(coords, contour),
     coords_(coords), val_(val), initialized_(true)
{
    // Delaunator is constructed in colon initialization
    if (val_.size() != coords_.size() /2){
        try {
            throw MeshException(1);
        } catch (MeshException& e) {
            std::cerr << "Found " << val.size() << " values for ";
            std::cerr << coords_.size() /2 << " pairs of coordinates." << std::endl;
            std::cout << e.what() << std::endl;
        }
    }
}

size_t TriagMesh::size()
{
    return val_.size();
}

size_t TriagMesh::numTriag()
{
    return d_.triangles.size()/3;
}

std::vector<size_t> TriagMesh::edgesOfTriag(size_t t)
{
    try {
        if (t > d_.triangles.size()){
            
            std::cout << t << " huh>>>" << std::endl;
            throw MeshException(4);
        }
    } catch (MeshException& e) {
        e.what();
    }
    t = t - (t % 3); // make sure we're at the beginning of the triangle
    std::vector<size_t> out {t, t +1, t + 2};
    return out;
}

std::vector<size_t> TriagMesh::pointsOfTriag(size_t t)
{
    std::vector<size_t> out(6);
    t = t - (t % 3); // make sure we're at the beginning of the triangle
    std::vector<size_t> edges = edgesOfTriag(t);
    for (int i = 0; i < 3; i++){
//        size_t coord_id = d_.triangles[t + i];
        size_t coord_id = d_.triangles.at(edges[i]); // safe retrieval
        out[2 * i]     = 2 * coord_id;     // x coordinate
        out[2 * i + 1] = 2 * coord_id + 1; // y coordinate
    }
    return out;
}

size_t TriagMesh::triagOfEdge(size_t e)
{
    try {
        if (e >= d_.triangles.size()){
            throw MeshException(5);
        }
    } catch (MeshException& e) {
        e.what();
    }
    return e - (e % 3);
}

size_t TriagMesh::nextHalfedge(size_t e)
{
    if (e % 3 == 2){
        return e - 2;
    } else {
        return e + 1;
    }
}

size_t TriagMesh::prevHalfedge(size_t e)
{
    if (e % 3 == 0){
        return e + 2;
    } else {
        return e - 1;
    }
}

std::vector<MeshPoint> TriagMesh::coordsOfTriag(size_t t)
{
    std::vector<MeshPoint> out(3);
    t = t - (t % 3); // go back to the beginning of the triangle
    std::vector<size_t> indices = pointsOfTriag(t);
    for (int i = 0; i < 3; i++){
        out[i].x_ = d_.coords.at(indices[ 2 * i]);
        out[i].y_ = d_.coords.at(indices[ 2 * i + 1]);
    }
    return out;
}

LineSeg TriagMesh::edgeToLineSeg(size_t e)
{
    size_t start = d_.triangles[e];
    double a_x = d_.coords.at(2 * start);
    double a_y = d_.coords.at(2 * start + 1);
    
    size_t end = d_.halfedges.at(e);
    size_t coord_id;
//    size_t end;
    if (end == delaunator::INVALID_INDEX){ // there's no opposite triangle
        // look for the start of next half edge
        coord_id = d_.triangles.at((e % 3 == 2) ? e - 2 : e + 1);
    }
    else {
        coord_id = d_.triangles.at(end);
    }
    double b_x = d_.coords.at(2 * coord_id);
    double b_y = d_.coords.at(2 * coord_id + 1);
    
    MeshPoint pa(a_x, a_y);
    MeshPoint pb(b_x, b_y);
    LineSeg rtn(pa, pb);
    return rtn;
}

size_t TriagMesh::neighborTriag(size_t e)
{
    size_t opposite = d_.halfedges.at(e);
    size_t triag;
    if (opposite != delaunator::INVALID_INDEX){
        triag = triagOfEdge(opposite);
    } else {
        triag = delaunator::INVALID_INDEX; // neighboring triangle doesn't exist
    }
    return triag;
}


VecDoub TriagMesh::barycentric(MeshPoint point, size_t t)
{
//    assert(point.size() == 2); // make sure that point is a valid coord
    double x = point.x_;
    double y = point.y_;
    
    std::vector<MeshPoint> coords = coordsOfTriag(t);
    double x1 = coords[0].x_;
    double y1 = coords[0].y_;
    double x2 = coords[1].x_;
    double y2 = coords[1].y_;
    double x3 = coords[2].x_;
    double y3 = coords[2].y_;
    
    double inv_det = 1/ ( (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3) );
    VecDoub out(3);
    out[0] = inv_det * ( (y2 - y3) * (x - x3) + (x3 - x2) * (y - y3) );
    out[1] = inv_det * ( (y3 - y1) * (x - x3) + (x1 - x3) * (y - y3) );
    out[2] = 1 - out[0] - out[1];
    
    return out;
}


bool TriagMesh::isInTriag(MeshPoint point, size_t t)
{
    VecDoub bar = barycentric(point, t);
    bool rtn = true;
    for (int i = 0; i < bar.size(); i++){
        rtn = rtn && (bar[i] >= 0) && (bar[i] <= 1);
    }
    return rtn;
}

MeshPoint TriagMesh::centroid(size_t t)
{
    std::vector<MeshPoint> coords = coordsOfTriag(t);
    double sumx(0), sumy(0);
    for (int i = 0; i < coords.size(); i++){
        sumx += coords[i].x_;
        sumy += coords[i].y_;
    }
    double xout = sumx/3;
    double yout = sumy/3;
    MeshPoint rtn(xout, yout);
    return rtn;
}


MeshPoint TriagMesh::circumcenter(size_t t)
{
    std::vector<MeshPoint> coords = coordsOfTriag(t);
    MeshPoint rtn;
    try {
        MeshPoint pa, pb, pc;
        pa = coords[0];
        pb = coords[1];
        pc = coords[2];
        std::pair<double, double> out = delaunator::circumcenter(pa.x_, pa.y_, pb.x_,  pb.y_, pc.x_, pc.y_);
        rtn.x_ = std::get<0>(out);
        rtn.y_ = std::get<1>(out);
    } catch (std::underflow_error & e) {
        // a silent nan if circumcenter underflew
        rtn.x_ = std::nan(0);
        rtn.y_ = std::nan(0);
    }
    return rtn;
}

void TriagMesh::edgesAroundPoint(size_t start, std::vector<size_t> &eout)
{
    // start is the starting t index
    size_t outgoing = d_.triangles[start]; // outgoing half edge at start
//    size_t origin = d_.triangles[d_.halfedges[start]]; // points toward center
    size_t origin = d_.halfedges[start]; // points toward center
    size_t incoming = origin;
    
//    size_t incoming = start;
    while (incoming != delaunator::INVALID_INDEX) {
        eout.push_back(incoming);
        size_t outgoing = nextHalfedge(incoming);
        incoming = d_.halfedges.at(outgoing);
        if (incoming == delaunator::INVALID_INDEX){
            throw MeshException(10);
        }
        if( incoming == origin) break;
    }
    
    
//    bool boundary = false;
//    // clockwise search
//    while (incoming != start){
//        eout.push_back(incoming);
//        size_t outgoing = nextHalfedge(incoming);
//        // check if we've hit the convex hull. Can't recover from that
//        if (outgoing == delaunator::INVALID_INDEX){
//            boundary  = true; // signal for hitting the boundary
//            break;
//        } else {
//            incoming = d_.halfedges.at(outgoing);
//        }
//    }
//    // counterclockwise search
//    if (boundary){
//        size_t newstart = d_.halfedges.at(start);
//        size_t outgoing = newstart;
//        while (outgoing != newstart){
//            incoming = prevHalfedge(outgoing);
//            if (out.find(incoming) != out.end()){
//                out.insert(incoming);
//            }
//            outgoing = d_.halfedges.at(incoming);
//            if (outgoing == delaunator::INVALID_INDEX){
//                break;
//            }
//        }
//    }
}

void TriagMesh::findVeroni(size_t start, std::vector<size_t> &out)
{
    std::vector<size_t> edges;
//    std::set<size_t> visited;
    try {
        edgesAroundPoint(start, edges); // all the edges around this point
        // find the triangles that include all these edges
        for (int i = 0; i < edges.size(); i++){
            size_t t = triagOfEdge(edges.at(i));
//            std::cerr << "hi! " << t << std::endl;
            out.push_back(t);
//            if (out.find(t) == out.end()){
////                std::cerr << "inserting! " << t << std::endl;
//                out.insert(t);
//            }
        }
    } catch (MeshException& e) {
        e.what();
    }
}

std::vector<size_t> TriagMesh::search(MeshPoint p, size_t& init)
{
    if (isRaw()){
        throw MeshException(0);
    }
    std::vector<size_t> rtn;
    size_t t_now(init);
    rtn.push_back(t_now);
//    size_t t_next;
    while (!isInTriag(p, t_now)){
//        std::cout << "searching t: " << t_now << std::endl;
        size_t t_next = walk(p, t_now);
        t_now = t_next;
        rtn.push_back(t_now);
    }
    init = t_now;
    return rtn;
}

size_t TriagMesh::walk(MeshPoint p, size_t t_now)
{
    std::vector<size_t> edges_head = edgesOfTriag(t_now);
    
    size_t intersection(delaunator::INVALID_INDEX); // an invalid result to initialize
    for (int i=0; i< edges_head.size(); i++){
        size_t e = edges_head[i];
        LineSeg edge = edgeToLineSeg(e);
        MeshPoint center = centroid(t_now);
        LineSeg query(p, center);
        if (query.isCross(edge)){
            intersection = e;
        }
    }
    if (intersection == delaunator::INVALID_INDEX){
        // no intersection found
        throw MeshException(2);
    }
    size_t e_opposite = d_.halfedges[intersection];
    if (e_opposite == delaunator::INVALID_INDEX){
        // point is outside domain
        throw MeshException(3);
    }
    size_t t_next = triagOfEdge(e_opposite);
    // program success
    return t_next;
    // function should terminate before it gets here
//    return delaunator::INVALID_INDEX;
}

double TriagMesh::interp(MeshPoint p, size_t& init)
{
    return interp(p, val_, init);
}

double TriagMesh::interp(MeshPoint p, VecDoub& val, size_t& init)
{
    if (val.size() != coords_.size() / 2){
        try {
            throw MeshException(1);
        } catch (MeshException& e) {
            std::cerr << "Found " << val.size() << " values for ";
            std::cerr << coords_.size() /2 << " pairs of coordinates." << std::endl;
            std::cout << e.what() << std::endl;
        }
    }
    std::vector<size_t> triag = search(p, init);
    VecDoub bary = barycentric(p, triag.back());
    
    std::vector<size_t> points = pointsOfTriag(triag.back());
    double out(0);
    for (int i=0; i<6; i+=2){
        size_t valIndex = points.at(i) / 2;
        double k = bary.at(i/2);
        double value = val.at(valIndex);
        out += k * value;
    }
    return out;
}

void TriagMesh::print(const char* fname){
    // print triangulation to file
    FILE * pFile;
    pFile = fopen (fname,"w");
    for(std::size_t i = 0; i < d_.triangles.size(); i+=3) {
        fprintf(pFile,
            "%f, %f\n %f, %f\n %f, %f\n",
            d_.coords[2 * d_.triangles[i    ]    ],    //tx0
            d_.coords[2 * d_.triangles[i    ] + 1],    //ty0
            d_.coords[2 * d_.triangles[i + 1]    ],    //tx1
            d_.coords[2 * d_.triangles[i + 1] + 1],    //ty1
            d_.coords[2 * d_.triangles[i + 2]    ],    //tx2
            d_.coords[2 * d_.triangles[i + 2] + 1]     //ty2
        );
    }
    fclose(pFile);
}

void TriagMesh::print(std::string& fname)
{
    // print triangulation to file
    std::ofstream pFile;
    pFile.open(fname);
    for(std::size_t i = 0; i < d_.triangles.size(); i+=3) {
        for (int j = 0; j < 3; ++j){
            size_t ix = 2 * d_.triangles[i + j];
            size_t iy = 2 * d_.triangles[i + j] + 1;
            pFile << d_.coords[ix] << ',';
            pFile << d_.coords[iy] << std::endl;;
        }
    }
    pFile.close();
}

void TriagMesh::print(std::string& fname, std::vector<int> datamask)
{
    // print triangulation to file
    assert(datamask.size() == d_.coords.size()/2);
    std::ofstream pFile;
    pFile.open(fname);
    for(std::size_t i = 0; i < d_.triangles.size(); i+=3) {
        size_t tx0 = d_.triangles[i    ];
        size_t tx1 = d_.triangles[i + 1];
        size_t tx2 = d_.triangles[i + 2];
//        int sum = datamask.at(tx0) + datamask.at(tx1) + datamask.at(tx2);
        
        int mask0 = datamask.at(tx0);
        int mask1 = datamask.at(tx1);
        int mask2 = datamask.at(tx2);
        // see if all 3 points are wall
        bool allwall = (mask0 == 2 && mask1 == 2 && mask2 == 2);
//        bool allcore = (mask0 == 0 && mask1 == 0 && mask2 == 0);
        if (!allwall){
            for (int j = 0; j < 3; ++j){
                size_t ix = 2 * d_.triangles[i + j];
                size_t iy = 2 * d_.triangles[i + j] + 1;
                pFile << d_.coords[ix] << ',';
                pFile << d_.coords[iy] << std::endl;;
            }
        }
    }
    pFile.close();
}

void TriagMesh::print(std::string& fname, VecDoub& data, std::vector<int> datamask)
{
    // print triangulation to file
    assert(datamask.size() == d_.coords.size()/2);
    std::ofstream pFile;
    pFile.open(fname);
    for(std::size_t i = 0; i < d_.triangles.size(); i+=3) {
        size_t tx0 = d_.triangles[i    ];
        size_t tx1 = d_.triangles[i + 1];
        size_t tx2 = d_.triangles[i + 2];
//        int sum = datamask.at(tx0) + datamask.at(tx1) + datamask.at(tx2);
        
        int mask0 = datamask.at(tx0);
        int mask1 = datamask.at(tx1);
        int mask2 = datamask.at(tx2);
        // see if all 3 points are wall
        bool allwall = (mask0 == 2 && mask1 == 2 && mask2 == 2);
//        bool allcore = (mask0 == 0 && mask1 == 0 && mask2 == 0);
        if (!allwall){
            for (int j = 0; j < 3; ++j){
                size_t ix = 2 * d_.triangles[i + j];
                size_t iy = 2 * d_.triangles[i + j] + 1;
                pFile << d_.coords[ix] << ',';
                pFile << d_.coords[iy] << ',';
                pFile << data.at(d_.triangles[i + j]) << std::endl;
            }
        }
    }
    pFile.close();
}

void TriagMesh::printVal(std::string& fname)
{
    // print mesh grid nodes and values to file
    std::ofstream pFile;
    pFile.open(fname);
    for(int i = 0; i < coords_.size(); i += 2){
        pFile << coords_.at(i) << ',' << coords_.at(i+1) << ',';
        pFile << val_.at(int(i/2)) << std::endl;
    }
    pFile.close();
}

void TriagMesh::printVal(const char* fname, VecDoub& data)
{
    std::string file(fname);
    printVal(file, data);
}

void TriagMesh::printVal(std::string& fname, VecDoub& data)
{
    // print mesh grid nodes and provided values to file
    std::ofstream pFile;
    pFile.open(fname);
    pFile << "xcoord, ycoord, value" << std::endl;
    for(int i = 0; i < coords_.size(); i += 2){
        pFile << coords_.at(i) << ',' << coords_.at(i+1) << ',';
        pFile << data.at(int(i/2)) << std::endl;
    }
    pFile.close();
}

bool TriagMesh::isRaw()
{
    return !initialized_;
}

void TriagMesh::printVeroni(std::string& prefix, std::vector<int> datamask)
{
    assert(datamask.size() == d_.coords.size()/2);
    std::ofstream p1, p2;
    p1.open(prefix + "veroniCoords.txt");
    p2.open(prefix + "veroniEdges.txt");

    // first print all the circumcenters to veroniNodes
//    p1 << "v_index, t_index, x, y" << std::endl;
    p1 << "//v_index, x, y" << std::endl;
    for (std::size_t v = 0; v < d_.triangles.size()/3; v++){
        std::size_t t = v * 3;
//        p1 << v << ',' << t << ',';
        p1 << v << ',';
        MeshPoint c = circumcenter(t);
        p1 << c.x_ << ',' << c.y_ << std::endl;
    }
    
    // now print all the connections
//    p2 << "t_index, v_index[multiple]" << std::endl;
    p2 << "//val_index, n_vertices, v_index[multiple]" << std::endl;
    std::set<size_t> visited;
    
    for (std::size_t i = 0; i < d_.triangles.size(); i++){
        size_t valIndex = d_.triangles[i];
        std::cerr << "valindex " << valIndex << std::endl;
        int mask = datamask.at(valIndex);
        if (mask != 2){ // if point is not on the wall
            std::vector<size_t> triags;
            
            size_t icoord_x = 2 * d_.triangles[i]; // this is not working
            std::cerr << "icoord: " << icoord_x << std::endl;
            
            if (visited.find(icoord_x) == visited.end()){ // we've never been to this point before
                std::cerr << "hi!" << std::endl;
                visited.insert(icoord_x); // mark that we've been here
                
                // now search for all the edges
                findVeroni(i, triags);
                // print to file
//                p2 << i << ',' << valIndex << ',';
                p2 << valIndex << ',' << triags.size() << ',';
                // valIndex will point us to the node in the value list
                std::vector<size_t>::iterator it;
                for (it = triags.begin(); it != triags.end(); ++it){
                    std::cerr << "found " << (*it)<< std::endl;
                    size_t v = (*it)/3;
                    p2 << v << ',';
                }
                p2 << std::endl;
            }
        }
    }
    p1.close();
    p2.close();
    return;
}

void TriagMesh::printTriag(std::string& prefix, std::vector<int>& datamask)
{
    std::ofstream p1, p2;
    p1.open(prefix + "triag.node");
    p2.open(prefix + "triag.ele");
    
    // first write .node files:
    // <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
    size_t n = coords_.size() / 2;
    p1 << "//n_vertices, n_dim, n_att, n_bndry" << std::endl;
    p1 << "//i_vertices, x, y, <att>, bndry" << std::endl;
    p1 << n << ",2,1,1" << std::endl;
    
    for (int i = 0; i < val_.size(); i++){
        bool bdry = (datamask.at(i) == 2); // whether it's wall
        p1 << i << ',';
        size_t ix = i * 2;
        p1 << coords_.at(ix) << ',' << coords_.at(ix + 1) << ',';
        p1 << val_.at(i) << ',' << bdry << std::endl;
    }
    p1.close();
    
    p2 << "//<# of triangles> <nodes per triangle> <# of attributes>" << std::endl;
    p2 << "//<triangle #> <node> <node> <node> ... [attributes]" << std::endl;
    
    size_t numTri = d_.triangles.size()/3;
    p2 << numTri << ",3,0" << std::endl;
    for (int i = 0; i < numTri; i++){
        int i_tri = i * 3;
        p2 << i <<',';
        p2 << d_.triangles[i_tri] << ',';
        p2 << d_.triangles[i_tri + 1] << ',';
        p2 << d_.triangles[i_tri + 2] << std::endl; // no surface quantities
    }
}

void TriagMesh::printTriag(std::string& prefix, std::vector<int>& datamask, std::vector<VecDoub>& values)
{
    std::ofstream p1, p2;
    p1.open(prefix + "triag.node");
    p2.open(prefix + "triag.ele");
    
    // first write .node files:
    // <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
    size_t n = coords_.size() / 2;
    p1 << "//n_vertices, n_dim, n_att, n_bndry" << std::endl;
    p1 << "//i_vertices, x, y, <att>, bndry" << std::endl;
    p1 << n << ",2,"<< values.size() <<",1" << std::endl;
    
    for (int i = 0; i < val_.size(); i++){
        bool bdry = (datamask.at(i) == 2); // whether it's wall
        p1 << i << ',';
        size_t ix = i * 2;
        p1 << coords_.at(ix) << ',' << coords_.at(ix + 1) << ',';
        
        for (int j = 0; j < values.size(); j++){
            p1 << values.at(j).at(i);
            if (j != values.size()){
                p1 << ',';
            }
        }
        p1 << bdry << std::endl;
    }
    p1.close();
    
    p2 << "//<# of triangles> <nodes per triangle> <# of attributes>" << std::endl;
    p2 << "//<triangle #> <node> <node> <node> ... [attributes]" << std::endl;
    
    size_t numTri = d_.triangles.size()/3;
    p2 << numTri << ",3,"<<values.size() << std::endl;
    for (int i = 0; i < numTri; i++){
        int i_tri = i * 3;
        p2 << i <<',';
        p2 << d_.triangles[i_tri] << ',';
        p2 << d_.triangles[i_tri + 1] << ',';
        p2 << d_.triangles[i_tri + 2] << ','; // no surface quantities
        
        size_t icoord = d_.triangles[i_tri];
        size_t icoord1 = d_.triangles[i_tri + 1];
        size_t icoord2 = d_.triangles[i_tri + 2];
//        size_t ival = icoord/2;
        
        for (int j = 0; j < values.size(); j++){
            double sum = values.at(j).at(icoord) + values.at(j).at(icoord1) + values.at(j).at(icoord2);
            double val = sum / 3.0;
            p2 << val;
            if (j != values.size() - 1){
                p2 << ',';
            } else {
                p2 << std::endl;
            }
        }
        
    }
}

void TriagMesh::printTriag(std::string& prefix, std::vector<int>& datamask, std::vector<bool>& valid, std::vector<VecDoub>& values)
{
    std::ofstream p1, p2;
    p1.open(prefix + "triag.node");
    p2.open(prefix + "triag.ele");
    
    // first write .node files:
    // <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
    size_t n = coords_.size() / 2;
    p1 << "//n_vertices, n_dim, n_att, n_bndry" << std::endl;
    p1 << "//i_vertices, x, y, <att>, bndry" << std::endl;
    p1 << n << ",2,"<< values.size() <<",1" << std::endl;
    
    for (int i = 0; i < val_.size(); i++){
        bool bdry = (datamask.at(i) == 2); // whether it's wall
        p1 << i << ',';
        size_t ix = i * 2;
        p1 << coords_.at(ix) << ',' << coords_.at(ix + 1) << ',';
        
        for (int j = 0; j < values.size(); j++){
            p1 << values.at(j).at(i);
            if (j != values.size()){
                p1 << ',';
            }
        }
        p1 << bdry << std::endl;
    }
    p1.close();
    
    // now write .ele files
    p2 << "//<# of triangles> <nodes per triangle> <# of attributes>" << std::endl;
    p2 << "//<triangle #> <node> <node> <node> ... [attributes]" << std::endl;
    
    size_t numTri = d_.triangles.size()/3;
    p2 << numTri << ",3," << values.size()+1 << std::endl;
    for (int i = 0; i < numTri; i++){
        int i_tri = i * 3;
        p2 << i <<',';
        
        size_t icoord0 = d_.triangles[i_tri];
        size_t icoord1 = d_.triangles[i_tri + 1];
        size_t icoord2 = d_.triangles[i_tri + 2];
        
        p2 << icoord0 << ',';
        p2 << icoord1 << ',';
        p2 << icoord2 << ','; // no surface quantities
        
//        size_t ival0 = icoord0/2;
//        size_t ival1 = icoord1/2;
//        size_t ival2 = icoord2/2;
        
        for (int j = 0; j < values.size(); j++){
            double val0, val1, val2;
//            val0 = values.at(j).at(ival0);
//            val1 = values.at(j).at(ival1);
//            val2 = values.at(j).at(ival2);
            val0 = values.at(j).at(icoord0);
            val1 = values.at(j).at(icoord1);
            val2 = values.at(j).at(icoord2);
            double sum = val0 + val1 + val2;
            double val = sum / 3.0;
            p2 << val << ',';
//            if (j != values.size() - 1){
//                p2 << ',';
//            } else {
//                p2 << std::endl;
//            }
        }
        p2 << valid.at(i_tri) << std::endl;
        
    }
}
