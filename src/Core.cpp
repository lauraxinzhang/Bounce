/**
 * \file    Core.cpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    October, 2020
 *
 * \brief   Implements the Core class
 * \details Stores the core plasma profile
 */

#include "Core.hpp"

Core::Core(std::string dirIn, const std::string& eqdsk)
//	 :trRMaj_(0), trTe_(0), trTi_(0), trNe_(0), trNi_(0)
     : dirIn_(dirIn), mTe_(nullptr), mTi_(nullptr), mNe_(nullptr), mNi_(nullptr), machine_(nullptr)
{
//    std::string gpath = dirIn + "g204202.00461";
    std::string gpath = dirIn + eqdsk;
    machine_ = new GEquilibrium(gpath);
    GEquilibrium nstx(gpath);
    VecDoub gPsiGrid(nstx.rGrid_.size());
    
    for (int i = 0; i < nstx.rGrid_.size(); ++i){
        MeshPoint p(nstx.rGrid_.at(i), 0);
        double psi = nstx.getPsiNorm(p);
        gPsiGrid.at(i) = psi;
    }
    XYMesh psiOfR(nstx.rGrid_, gPsiGrid);
    
    VecDoub Te = readTrData("TE");
    VecDoub Ti = readTrData("TI");
    VecDoub Ne = readTrData("NE");
    VecDoub Ni = readTrData("ND");
    VecDoub RMaj = readTrData("RMAJM");
    
    VecDoub trPsi(Te.size()); // value of normalized POLOIDAL flux
    for (int i = 0; i < Te.size(); i++){
        double Rcenter = 0.5 * ( RMaj.at(Te.size() + i) + RMaj.at(Te.size() + i + 1) );
        MeshPoint p(Rcenter/100, 0); // convert to meter
        //TODO: add function in Mesh to handle single value input, instead of MeshPoint
        double psi = psiOfR.interp(p);
        trPsi.at(i) = psi;
        
        // convert unites for density
        Ne.at(i) = Ne.at(i) * 1E6;
        Ni.at(i) = Ni.at(i) * 1E6;
    }
    
    trTe_ = Te;
    trTi_ = Ti;
    trNe_ = Ne;
    trNi_ = Ni;
    trRMaj_ = RMaj;
    trPsi_ = trPsi;
    
    // plasma parameters interpolated on normalized POLOIDAL flux
    mTe_ = new XYMesh(trPsi, Te);
    mTi_ = new XYMesh(trPsi, Ti);
    mNe_ = new XYMesh(trPsi, Ne);
    mNi_ = new XYMesh(trPsi, Ni);
    
}

void Core::setDirIn(std::string& dirIn)
{
    dirIn_ = dirIn;
    return;
}

VecDoub Core::readTrData(const char *  var)
{
    Parser parse;
//    std::string trDir = "204202R90_461/";
    std::string suffix = ".txt";
    std::string fname = var + suffix;
    VecDoub out = parse.parseVec(dirIn_, fname);
    return out;
}

double Core::getTemp(double psiPol, Particle background)
{
    double rtn(-1); // invalid value
    MeshPoint p(psiPol, 0);
    try {
        if (background.charge() < 0){ // electron
            if (psiPol > 0 && psiPol <= trPsi_.front()){
                rtn = trTe_.front();
            } else if(psiPol >= trPsi_.back() && psiPol <= 1) {
                rtn = trTe_.back();
            } else {
                rtn =  mTe_->interp(p);
            }
        } else {
            if (psiPol > 0 && psiPol <= trPsi_.front()){
                rtn = trTi_.front();
            } else if(psiPol >= trPsi_.back() && psiPol <= 1) {
                rtn = trTi_.back();
            } else {
                rtn =  mTi_->interp(p);
            }
        }
    } catch (MeshException& e) {
        if (e.code_ == 9) {
            // query point maybe in SOL, return a silent nan
            rtn = std::nan("0");
        } else {
            e.what();
        }
    }
    return rtn;
}

double Core::getDense(double psiPol, Particle background)
{
    double rtn(-1); // invalid value
    try {
        MeshPoint p(psiPol, 0);
        if (background.charge() < 0){ // electron
            if (psiPol > 0 && psiPol <= trPsi_.front()){
                rtn = trNe_.front();
            } else if(psiPol >= trPsi_.back() && psiPol <= 1) {
                rtn = trNe_.back();
            } else {
                rtn =  mNe_->interp(p);
            }
        } else {
            if (psiPol > 0 && psiPol <= trPsi_.front()){
                rtn = trNi_.front();
            } else if(psiPol >= trPsi_.back() && psiPol <= 1) {
                rtn = trNi_.back();
            } else {
                rtn =  mNi_->interp(p);
            }
        }
    } catch (MeshException& e) {
        if (e.code_ == 9) {
            // query point maybe in SOL, return a silent nan
            rtn = std::nan("0");
        } else {
            e.what();
        }
    }
    return rtn;
}

double Core::getTemp(Vector pos, Particle background, size_t& t_init)
{
    double rr = sqrt(pos.x() * pos.x() + pos.y() * pos.y() );
    double zz = pos.z();
    
    MeshPoint p(rr, zz);
    // cast to Psi
    double psiNorm = machine_ -> getPsiNorm(p);
    return getTemp(psiNorm, background);
}

double Core::getDense(Vector pos, Particle background, size_t& t_init)
{
    double rr = sqrt(pos.x() * pos.x() + pos.y() * pos.y() );
    double zz = pos.z();
    
    MeshPoint p(rr, zz);
    // cast to Psi
    double psiNorm = machine_ -> getPsiNorm(p);
    return getDense(psiNorm, background);
}

void Core::getValSep(double& Te, double& Ti, double& n0)
{
    Te = trTe_.back();
    Ti = trTi_.back();
    n0 = trNe_.back();
    return;
}
