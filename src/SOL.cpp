/**
 * \file    SOL.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2020
 *
 * \brief   Implements the SOL class
 */

#include "SOL.hpp"

//SOL::SOL()
//    :dirOut_("/Users/xinzhang/Documents/thesis/SOLFI/output/")
//{
////    std::cout << "isRaw? " << mesh_.isRaw() << std::endl;
//}
//
//SOL::SOL(VecDoub& coords, std::string dirOut)
//{
////    std::cout << mesh_.isRaw() << std::endl;
//}
//
//void SOL::setDirOut(std::string& dirOut)
//{
//    dirOut_ = dirOut;
//}
//
//void SOL::printMeshTriag()
//{
//    std::string fname = dirOut_ + "SOLtriag.txt";
////    mesh_->print(fname);
//    std::cerr << "TODO: Printing in SOL not implemented" << std::endl;
//}

SOLPS::SOLPS(const std::string& dirIn)
    :dirOut_("./output/")
{
    Parser parser;
    std::string solpsIn = dirIn + "solpsprofiles/";
    MatDoub crx0 = parser.parseDatMat(solpsIn , "crx0.dat");
    MatDoub crx1 = parser.parseDatMat(solpsIn , "crx1.dat");
    MatDoub crx2 = parser.parseDatMat(solpsIn , "crx2.dat");
    MatDoub crx3 = parser.parseDatMat(solpsIn , "crx3.dat");
    std::vector<MatDoub> crx = {crx0, crx1, crx2, crx3};

    MatDoub cry0 = parser.parseDatMat(solpsIn , "cry0.dat");
    MatDoub cry1 = parser.parseDatMat(solpsIn , "cry1.dat");
    MatDoub cry2 = parser.parseDatMat(solpsIn , "cry2.dat");
    MatDoub cry3 = parser.parseDatMat(solpsIn , "cry3.dat");
    std::vector<MatDoub> cry = {cry0, cry1, cry2, cry3};
    
    MatDoub lix = parser.parseDatMat(solpsIn , "leftix.dat");
    MatDoub rix = parser.parseDatMat(solpsIn , "rightix.dat");
    MatDoub tix = parser.parseDatMat(solpsIn , "topix.dat");
    MatDoub bix = parser.parseDatMat(solpsIn , "bottomix.dat");
    std::vector<MatDoub> nbx = {lix, tix, rix, bix};
    
    MatDoub liy = parser.parseDatMat(solpsIn , "leftiy.dat");
    MatDoub riy = parser.parseDatMat(solpsIn , "rightiy.dat");
    MatDoub tiy = parser.parseDatMat(solpsIn , "topiy.dat");
    MatDoub biy = parser.parseDatMat(solpsIn , "bottomiy.dat");
    std::vector<MatDoub> nby = {liy, tiy, riy, biy};
    
    MatDoub ne = parser.parseDatMat(solpsIn , "ne.dat");
    MatDoub te = parser.parseDatMat(solpsIn , "te_eV.dat");
    MatDoub nD = parser.parseDatMat(solpsIn , "nD1.dat");
    MatDoub ti = parser.parseDatMat(solpsIn , "ti_eV.dat");
    std::vector<MatDoub> vals = {ne, te, nD, ti};
    
    std::cout << ne.getW() << ',' << ne.getH() << std::endl;
    
    setupTriag(crx, cry);
    quad_ = new QuadMesh(crx, cry, nbx, nby, vals);
}

void SOLPS::setupTriag(std::vector<MatDoub>& crx, std::vector<MatDoub>& cry)
{
    nx_ = crx.front().getW();
    ny_ = cry.front().getH();
    VecDoub * coords = new VecDoub(nx_ * ny_ * 2);
    VecDoub val(nx_ * ny_, 0);
    for (int iy = 0; iy < ny_; ++iy) {
        for (int ix = 0; ix < nx_; ++ix){
            // lower left corner
            MeshPoint ll(crx.at(0).at(ix, iy), cry.at(0).at(ix, iy));
            // lower right corner
            MeshPoint lr(crx.at(1).at(ix, iy), cry.at(1).at(ix, iy));
            // upper left corner
            MeshPoint ul(crx.at(2).at(ix, iy), cry.at(2).at(ix, iy));
            // upper right corner
            MeshPoint ur(crx.at(3).at(ix, iy), cry.at(3).at(ix, iy));

            Quad q(ll, lr, ul, ur);
            MeshPoint centroid = q.centroid();

            int i = iy * nx_ + ix; // flatten the matrix
            coords->at(2 * i)     = centroid.x_;
            coords->at(2 * i + 1) = centroid.y_;
        }
    }
    triag_ = new TriagMesh(*coords, val);
    return;
}

int SOLPS::getW()
{
    return nx_;
}

int SOLPS::getH()
{
    return ny_;
}

VecDoub SOLPS::getVal(MeshPoint& p, size_t& t_cache, bool histOut)
{
    // starting from the current cached t as initial guess
    std::vector<size_t> triagHist = triag_ -> search(p, t_cache);
    // final result of the triangulation search
    size_t i_triag = triagHist.back();
//    std::cerr << "found in triangle: " << i_triag << std::endl;
     // modify the input parameter t_cache to temporily store the
     // result of this current search, for use in the next call.
    t_cache = i_triag;
    std::vector<size_t> points = triag_ -> pointsOfTriag(i_triag);
    
    // some initialization for search on quad
    VecDoub vals;
    std::vector<size_t> q_hist;
    bool success(false);
    for (int i_flat = 0; i_flat < 3; ++i_flat){
        size_t vertex = points.at(2 * i_flat) / 2;
        size_t ix = (int) vertex % nx_; // remainder
        size_t iy = (int) vertex / nx_; // integer division
        std::vector<size_t> q_init = {ix, iy};
        try {
//            std::cout << "i_flat: " << i_flat << std::endl;
//            std::cout << "vertex: " << vertex << std::endl;
//            std::cout << "ix, iy: " << ix <<',' << iy << std::endl;
            q_hist = quad_ -> search(p, q_init);
            vals = quad_ -> interp(p, q_init);
//            std::cout << "value found in quad:" << std::endl;
            Quad q = quad_->ConstructQuad(q_init);
//            std::cout << q;
            success = true;
            break;
        } catch (MeshException& e) {
            // It is possible to recover if one of the triangle vertex leads to a deadend.
            assert(e.code_ == 7);
//            e.what();
            continue;
        }
    }
    try {
        if (!success){
            throw SOLException(2);
        } else {
            if (histOut) {
                printHist(triagHist, q_hist);
            }
        }
    } catch (SOLException& e) {
//        e.what();
    }
    return vals;
}

double SOLPS::getTemp(Vector pos, Particle background, size_t& t_init)
{
    double rr = sqrt(pos.x() * pos.x() + pos.y() * pos.y() );
    double zz = pos.z();
    MeshPoint p(rr, zz);
    VecDoub out = getVal(p, t_init);
    if (background.charge() < 0) {
        return out.at(1); // electron temperature
    } else {
        return out.at(3) * 0.75; // ion temperature artificially shifted
    }
}

double SOLPS::getDense(Vector pos, Particle background, size_t& t_init)
{
    double rr = sqrt(pos.x() * pos.x() + pos.y() * pos.y() );
    double zz = pos.z();
    MeshPoint p(rr, zz);
    VecDoub out = getVal(p, t_init);
    try {
        if (background.charge() < 0) {
            return out.at(0); // electron density
        } else {
            return out.at(2);
        }
    } catch (std::out_of_range& oor) { // is out is empty
        throw MeshException(3);
    }
//    if (spec) {
//        return out.at(0); // electron density
//    } else {
//        return out.at(2);
//    }
}

void SOLPS::printHist(std::vector<size_t> &t_hist, std::vector<size_t> &q_hist)
{
    std::string tfname = dirOut_ + "hist_t.txt";
    std::string qfname = dirOut_ + "hist_q.txt";
    std::string qcfname = dirOut_ + "coords_q.txt";
    std::ofstream tfile, qfile, qcfile;
    
    try {
        tfile.open(tfname, std::ios_base::app);
        qfile.open(qfname, std::ios_base::app);
        qcfile.open(qcfname, std::ios_base::app);
        if (!tfile.is_open() || !qfile.is_open() || !qcfile.is_open()){
            throw SOLException(3);
        } else {
//            for(int i=0; i < t_hist.size(); ++i){
//                tfile << t_hist.at(i) << std::endl;
//            }
//            for(int i=0; i < q_hist.size(); ++i){
//                size_t ix, iy;
//                ix = q_hist.at(i);
//                iy = q_hist.at(i+1);
//                qfile << ix << ',' << iy << std::endl;
//
//                //print the coordinates of the quads
//                std::vector<size_t> q_index = {ix, iy};
//                Quad q = quad_ -> ConstructQuad(q_index);
//                qcfile << q;
//                ++i;
//            }
            tfile << t_hist.size() << std::endl;
            qfile << q_hist.size() / 2 << std::endl;
        }
    } catch (SOLException& e) {
        e.what();
    }
}

void SOLPS::printQuads()
{
    std::string fname = dirOut_ + "quads.txt";
    std::ofstream quads;
    try {
        quads.open(fname);
        if (!quads.is_open()){
            throw SOLException(3);
        } else {
            for (size_t iy = 0; iy < ny_; ++iy){
                for (size_t ix = 0; ix < nx_; ++ix){
                    std::vector<size_t> q_index = {ix, iy};
                    Quad q = quad_ -> ConstructQuad(q_index);
                    quads << q;
                }
            }
        }
    } catch (SOLException& e) {
        std::cerr << "in printQuads." << std::endl;
        e.what();
    }
}
