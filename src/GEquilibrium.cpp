/**
 * \file    GEquilibrium.cpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    December, 2020
 *
 * \brief   Implements the Tokamak Equilirbium class.
 */

#include "GEquilibrium.hpp"

GEquilibrium::GEquilibrium()
        : gfile_(0, 0), Br_(0, 0), Bz_(0, 0), Btor_(0, 0),
            rGrid_(0), zGrid_(0), siGrid_(0),
            mBr_(nullptr), mBz_(nullptr), mBtor_(nullptr), mPsiRZ_(nullptr)
{}

GEquilibrium::GEquilibrium(int nw, int nh)
        : gfile_(nw, nh), Br_(nw, nh), Bz_(nw, nh), Btor_(nw, nh),
            rGrid_(nw), zGrid_(nh), siGrid_(nw),
            mBr_(nullptr), mBz_(nullptr), mBtor_(nullptr), mPsiRZ_(nullptr)
{}

GEquilibrium::GEquilibrium(const std::string& path)
    :gfile_(0, 0), Br_(0, 0), Bz_(0, 0), Btor_(0, 0),
     rGrid_(0), zGrid_(0), siGrid_(0)
{
    Parser parse;
    Parser::eqdsk gfile = parse.parseEqdsk(path);
    gfile_ = gfile;
    rGrid_ = gfile.rGrid_;
    zGrid_ = gfile.zGrid_;
    siGrid_ = gfile.siGrid_;
    mPsiRZ_ = new UniMesh(rGrid_, zGrid_, gfile_.psirz_);
    maxis_ = MeshPoint(gfile.rmaxis, gfile.zmaxis);
    setupB();
    setupLimiter();
    return;
}

void GEquilibrium::setupB()
{
    int nw = gfile_.nw;
    int nh = gfile_.nh;
    Br_.resize(nw, nh);
    Bz_.resize(nw, nh);
    Btor_.resize(nw, nh);
//    std::cout << "size:" << Br_.getW() << std::endl;
    
    MatDoub Dr = Dij(gfile_.rGrid_);
    MatDoub Dz = Dij(gfile_.zGrid_);
    
    for (int i = 0; i < nw; i++){
        double Ri = gfile_.rGrid_.at(i);
//        double norm = -1 / (2 * PI * Ri);
        double norm = 1 / (Ri); // compliant with convention. this is correct
//        double norm = -2 * PI  / ( Ri);
//        double norm = 1;
        
        
        for (int j = 0; j < nh; j++){
            
            // calculate Br
            double rdotted(0);
            for (int l = 0; l < nh; l++){
                rdotted += Dz.at(l, j) * gfile_.psirz_.at(i, l);
            }
            Br_.at(i, j) = norm * rdotted;
            
            // calculate Bz
            double zdotted(0);
            for (int l = 0; l < nw; l++){
                zdotted += Dr.at(l, i) * gfile_.psirz_.at(l, j);
            }
            Bz_.at(i, j) = -1 * norm * zdotted;
            
            // calculate Btor
            double sumNumer(0);
            double sumDenom(0);
            double psiNow = gfile_.psirz_.at(i, j);
            int nsi = gfile_.siGrid_.size();
            
            for (int k = 0; k < nsi; k++){
                int gam = gamma(k);
                double psik = gfile_.siGrid_.at(k);
                double denom, numer;
                if ( fabs(psik - psiNow) < 1E-10) {
                    // We're trying to evaluate on grid points
                    Btor_.at(i, j) = gfile_.fpol_.at(k) / Ri;
                    break;
                } else {
                    denom = gam / (psiNow - psik);
                    numer = denom * gfile_.fpol_.at(k);
                    sumDenom += denom;
                    sumNumer += numer;
                }
            }
            Btor_.at(i, j) = sumNumer / (Ri * sumDenom);
        }
    }
    mBr_ = new UniMesh(rGrid_, zGrid_, Br_);
    mBz_ = new UniMesh(rGrid_, zGrid_, Bz_);
    mBtor_ = new UniMesh(rGrid_, zGrid_, Btor_);
    mRPsi_ = new UniMesh(rGrid_, siGrid_);
    return;
}

void GEquilibrium::setupLimiter()
{
    // setup limiter structure, find rlim on both sides
    Parser::limiter limit = gfile_.getLimiter();
    for (int i = 0; i < limit.limitr; i++){
        double rlim1 = limit.rlim_.at(i);
        double zlim1 = limit.zlim_.at(i);
        MeshPoint p1(rlim1, zlim1);
        MeshPoint p2;
        if (i != limit.limitr - 1){
            double rlim2 = limit.rlim_.at(i + 1);
            double zlim2 = limit.zlim_.at(i + 1);
            p2 = MeshPoint(rlim2, zlim2);
        } else { // if we're at the end, wrap around
            double rlim2 = limit.rlim_.at(0);
            double zlim2 = limit.zlim_.at(0);
            p2 = MeshPoint(rlim2, zlim2);
        }
        if (zlim1 * p2.y_ <= 0){// crossing z = 0
            // take average in case they're not the same
            double rAve = (rlim1 + p2.x_) / 2;
            if (rlim1 < gfile_.rmaxis){
                rlimHFS_ = rAve;
            } else {
                rlimLFS_ = rAve;
            }
        }
        LineSeg line(p1, p2);
        limiter_.push_back(line);
    }
    // setup separatrix, find rsep on both sides and (fake) xpoints
    double zmin(0), zmax(0);
    double rzmin(0),rzmax(0);
    for (int i = 0; i < limit.nbbbs; i++){
        double rbbbs1 = limit.rbbbs_.at(i);
        double zbbbs1 = limit.zbbbs_.at(i);
        if (zbbbs1 > zmax){
            zmax = zbbbs1;
            rzmax = rbbbs1;
        }
        if (zbbbs1 < zmin) {
            zmin = zbbbs1;
            rzmin = rbbbs1;
        }
        MeshPoint p1(rbbbs1, zbbbs1);
        MeshPoint p2;
        if (i != limit.nbbbs - 1){
            double rlim2 = limit.rbbbs_.at(i + 1);
            double zlim2 = limit.zbbbs_.at(i + 1);
            p2 = MeshPoint(rlim2, zlim2);
        } else { // if we're at the end, wrap around
            double rlim2 = limit.rbbbs_.at(0);
            double zlim2 = limit.zbbbs_.at(0);
            p2 = MeshPoint(rlim2, zlim2);
        }
        if (zbbbs1 * p2.y_ <= 0){// crossing z = 0
            // take average in case they're not the same
            double rAve = (rbbbs1 + p2.x_) / 2;
            if (rbbbs1 < gfile_.rmaxis){
                rsepHFS_ = rAve;
            } else {
                rsepLFS_ = rAve;
            }
        }
        LineSeg line(p1, p2);
        separatrix_.push_back(line);
    }
    xtop_ = MeshPoint(rzmax, zmax);
    xbot_ = MeshPoint(rzmin, zmin);
    return;
}

Vector GEquilibrium::getB(MeshPoint& p)
{
    // no initial guess is needed for uniform mesh
    double Br = mBr_ -> interp(p);
    double Bz = mBz_ -> interp(p);
    double Btor = mBtor_ -> interp(p);
    Vector rtn(Br, Btor, Bz);
    return rtn;
}

Vector GEquilibrium::getB(Vector& v)
{
    double rr = sqrt(v.x() * v.x() + v.y() * v.y());
    MeshPoint p(rr, v.z());
    Vector B_cyl = getB(p);
    Vector B_cart;
    B_cyl.cyl2Cart(v, B_cart);
    return B_cart;
}

double GEquilibrium::getPsi(MeshPoint& p)
{
    double rtn = mPsiRZ_ -> interp(p);
    return rtn;
}


double GEquilibrium::getPsiNorm(MeshPoint &p)
{
    double psi_raw = getPsi(p);
//    std::cerr << "siraw: " << psi_raw << std::endl;
//    std::cerr << "sibdry: " << gfile_.sibry << std::endl;
//    std::cerr << siGrid_.front() << ',' << siGrid_.back() << std::endl;
    double extent = fabs(gfile_.sibry - gfile_.simag);
    double rtn = sqrt(fabs(psi_raw - gfile_.simag)/extent);
    return rtn;
}

double GEquilibrium::getPsiNorm(Vector& v)
{
    double rr = sqrt(v.x() * v.x() + v.y() * v.y());
    MeshPoint p(rr, v.z());
    return getPsiNorm(p);
}

double GEquilibrium::getCosTheta(Vector& v)
{
    double rNow = sqrt( v.x() * v.x() + v.y() * v.y());
    double zNow = v.z();
    // Calculate the current cosine(theta), theta is poloidal angle
    double rdiff = rNow - gfile_.rmaxis;
    double zdiff = zNow - gfile_.zmaxis;
    double rMinor = sqrt( zdiff * zdiff + rdiff * rdiff);
    double costheta = rdiff / rMinor;
    return costheta;
}

bool GEquilibrium::isLimiter(Vector v)
{
    //std::cerr << "hello limiter" << std::endl;
    double rr = sqrt(v.x() * v.x() + v.y() * v.y());
    MeshPoint p(rr, v.z());
    MeshPoint p_prime(gfile_.rmaxis, gfile_.zmaxis);
    LineSeg point(p, p_prime);
    
    bool crosses(false);
    for (int i = 0; i < limiter_.size(); i++){
        LineSeg wall = limiter_.at(i);
        if (wall.isCross(point)){
            crosses = true;
            break;
        }
    }
    return crosses;
}

bool GEquilibrium::isLimiter(Vector v, LineSeg& wallSeg)
{
    //std::cerr << "hello limiter" << std::endl;
    double rr = sqrt(v.x() * v.x() + v.y() * v.y());
    MeshPoint p(rr, v.z());
    MeshPoint p_prime(gfile_.rmaxis, gfile_.zmaxis);
    LineSeg point(p, p_prime);
    
    bool crosses(false);
    for (int i = 0; i < limiter_.size(); i++){
        LineSeg wall = limiter_.at(i);
        if (wall.isCross(point)){
            crosses = true;
            wallSeg = wall;
            break;
        }
    }
    return crosses;
}

bool GEquilibrium::isSeparatrix(Vector v)
{
    //std::cerr << "hello limiter" << std::endl;
    double rr = sqrt(v.x() * v.x() + v.y() * v.y());
    MeshPoint p(rr, v.z());
    MeshPoint p_prime(gfile_.rmaxis, gfile_.zmaxis);
    LineSeg point(p, p_prime);
    
    bool crosses(false);
    for (int i = 0; i < separatrix_.size(); i++){
        LineSeg wall = separatrix_.at(i);
        if (wall.isCross(point)){
            crosses = true;
            break;
        }
    }
    return crosses;
}

bool GEquilibrium::isSeparatrix(MeshPoint p)
{
    MeshPoint p_prime(gfile_.rmaxis, gfile_.zmaxis);
    LineSeg point(p, p_prime);
    
    bool crosses(false);
    for (int i = 0; i < separatrix_.size(); i++){
        LineSeg wall = separatrix_.at(i);
        if (wall.isCross(point)){
            crosses = true;
            break;
        }
    }
    return crosses;
}

int GEquilibrium::gamma(int k)
{
    int rtn;
    if (k % 2 == 0) {
        rtn = 1;
    } else {
        rtn = -1;
    }
    return rtn;
}

MatDoub GEquilibrium::Dij(VecDoub& grid)
{
    double extent = grid.back() - grid.front();
    int nn = grid.size();

    // D is always a square matrix
    MatDoub D(nn, nn);
    
    // fill in off-diagnoal elements
    double gammaJ, gammaI;
    for (int i = 0; i < nn; i++){
        gammaI = gamma(i);
        if (i == 0 || i == nn-1 ){
            gammaI = gammaI * 0.5; // 2nd lowest order interpolation
        }
        for (int j = 0; j < nn; j++){
            gammaJ = gamma(j);
            if (j == 0 || j == nn-1 ){
                gammaJ = gammaJ * 0.5;
            }
            
            if (i != j) {
                double sign = gammaJ / gammaI;
                double val = (nn - 1) / (extent * (i-j) );
                D.at(j, i) = sign * val;
            } else {
                D.at(j, i) = 0;
            }
        }
    }
    // fill in on-diagnoal elements
    for (int i = 0; i < nn; i++){
        double sum(0);
        for (int m = 0; m < nn; m++){
            sum += D.at(m, i);
        }
        D.at(i, i) = -1 * sum;
    }
    return D;
}
