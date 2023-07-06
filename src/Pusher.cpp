/**
 * \file    Pusher.cpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    October, 2020
 *
 * \brief   Implements the Pusher class
 */

#include "Pusher.hpp"

Pusher::Pusher(Fields& fields)
        : plasma_(nullptr)
{
    fields_ = &fields;
}

Pusher::Pusher(Plasma& plasma, Fields& fields)
{
    plasma_ = &plasma;
    fields_ = &fields;
}

Pusher::Pusher(const std::string& dirIn, const std::string& eqdsk)
    : fields_(nullptr), plasma_(nullptr), partList_(0),
     dirOut_("./output/")
	{
        plasma_  = new Tokamak(dirIn, eqdsk);
        std::string gpath = dirIn + eqdsk;
        fields_ = new GEquilibrium(gpath);
//        sol_ = new SOLPS(dirIn);
    }

//Pusher::~Pusher()
//{
//    if (plasma_){
//        delete plasma_;
//    }
//    if (fields_) {
//        delete fields_;
//    }
//}

void Pusher::setDirOut(std::string out)
{
    dirOut_ = out;
}

void Pusher::initParts(int nparts, double energy, Vector pos, int mu, int Z)
{
    partList_.resize(nparts);
    Particle sample(pos, pos, mu, Z);
    double B0 = fields_->getB(pos).mod();
    double vt = sqrt(2 * energy * QE / sample.mass());
    
    std::default_random_engine generator(int(time(NULL)));
    std::normal_distribution<double> distribution(0, vt); // generate a Gaussian distributed velocity
    
    for (int i = 0; i < nparts; i++){
        double vx = distribution(generator); // generate 3 normal distributed velocities.
        double vy =  distribution(generator);
        double vz = distribution(generator);
        
        Vector vel(vx, vy, vz);
        Particle part(pos, vel, mu, Z);
        double dt = part.setDt(B0);
//        std::cerr << "dt: " << dt << std::endl;
        partList_.at(i) = part;
    }
    return;
}

void Pusher::pushAll()
{//TODO: not implemented
    return;
}

void Pusher::pushSingle(Particle& part, int iter, bool collide, std::string prefix, bool write)
{
    if (collide && !plasma_){
        std::cerr << "Warning: no plasma information found for thermal collisions." << std::endl;
    }
    std::ofstream coord = safeOpen(prefix + "_coords.txt");
    std::ofstream coordRZ = safeOpen(prefix + "_coords_RZ.txt");
    std::ofstream energy = safeOpen(prefix + "_energy.txt");
    std::ofstream time = safeOpen(prefix + "_time.txt");
    
    Vector pos = part.pos();
    double B0 = fields_->getB(pos).mod();
    double dt = part.setDt(B0);
    double tnow = 0;
    
    std::cerr << dt << std::endl;
    
    
    for (int i = 0; i < iter; i++){
        
        if (write && (i % 3 == 0)){
//                std::cerr << part.pos() << std::endl;
            coord << part.pos() << std::endl;
            double x = part.pos().x();
            double y = part.pos().y();
            double r = sqrt( x*x + y*y );
            coordRZ << r << ',' << part.pos().z() << std::endl;
            
            double enow = part.energy();
            energy << enow << std::endl;
            time << tnow << std::endl;
        }
        
        Vector posNow = part.pos();
        Vector BNow;
        Vector ENow(0, 0, 0); // no electric field
        try {
            BNow = fields_ -> getB(posNow);
        } catch (MeshException& e) {
            if (e.code_ == 9){
                std::cerr << "Can't find magnetic info" << std::endl;
                std::cerr << part.pos() << std::endl;
                part.lost();
                break;
            } else {
                e.what();
            }
        }
//        std::cerr << "01 " << part.pos() << std::endl;
//        std::cerr << "02 " << BNow << std::endl;
//        part.moveCyl(ENow, BNow);
        part.move(ENow, BNow);
//        std::cerr << "02 " << part.pos() << std::endl;
        tnow += dt;
        if (fields_ -> isLimiter(part.pos())){
            part.lost();
            std::cerr << "Particle lost to limiter" << std::endl;
            break;
        }
        
//        std::cerr << part.pos() << std::endl;
        if (collide && plasma_){
            double nu_s, root_D_para, root_D_perp;
            try {
//                std::cerr << "1 " << part.pos() << std::endl;
                plasma_ -> thermalCoef(part, nu_s, root_D_para, root_D_perp);// fills in the values of thermal coefficients
//                std::cerr << "2 " << part.pos() << std::endl;
            } catch (MeshException& e) {
//                std::cerr << "33 " << part.pos() << std::endl;
                if (e.code_ == 3 || e.code_ == 8 || e.code_ == 9){
//                    std::cerr << "can't find plasma to collide with. skipping collision" << std::endl;
//                    std::cerr << "3 " << part.pos() << std::endl;
//                    std::cerr << "3 vel " << part.vel() << std::endl;
                    continue;
                } else {
                    
                    e.what();
                }
            }
//            std::cerr << "4 " << part.pos() << std::endl;
//            std::cerr << "4 vel " << part.vel() << std::endl;
            part.collide(nu_s, root_D_para, root_D_perp);
//            std::cerr << "5 " << part.pos() << std::endl;
//            std::cerr << "5 vel " << part.vel() << std::endl;
        }
            
//        std::cerr << "hello" << std::endl;

            
        
    }
}

//void Pusher::thermalCoef(Particle part, double& nu_s, double& root_D_para, double& root_D_perp)
//{
//    Particle electron;
//    Particle ion;
//    ion.setSpec(0, mu_);
//
//    Vector pos = part.pos();
//    double psiNow = machine_ -> getPsiNorm(pos);
//
//    double xe, nu_e, xi, nu_i; // constants for electrons and main ions
//
//
//    if (psiNow <= INTERFACE){
//        core_ -> xAndNu_ab(part, electron, xe, nu_e);
//        core_ -> xAndNu_ab(part, ion, xi, nu_i);
//
//        nu_s        = core_ -> nu_s(part, electron, xe, nu_e) + \
//                      core_ -> nu_s(part, ion, xi, nu_i);
//        root_D_para = core_ -> root_D_para(part, xe, nu_e) + \
//                      core_ -> root_D_para(part, xi, nu_i);
//        root_D_perp = core_ -> root_D_perp(part, xe, nu_e) + \
//                      core_ -> root_D_perp(part, xi, nu_i);
//    } else {
//        sol_ -> xAndNu_ab(part, electron, xe, nu_e);
//        sol_ -> xAndNu_ab(part, ion, xi, nu_i);
//
//        nu_s        = sol_ -> nu_s(part, electron, xe, nu_e) + \
//                      sol_ -> nu_s(part, ion, xi, nu_i);
//        root_D_para = sol_ -> root_D_para(part, xe, nu_e) + \
//                      sol_ -> root_D_para(part, xi, nu_i);
//        root_D_perp = sol_ -> root_D_perp(part, xe, nu_e) + \
//                      sol_ -> root_D_perp(part, xi, nu_i);
//    }
//    return;
//}

std::ofstream Pusher::safeOpen(std::string fname)
{
    std::string path = dirOut_ + fname;
    std::ofstream myfile(path); // open eqdsk file

    if (!myfile) {
        std::cerr << "Pusher: Unable to open file: " + fname << std::endl;
        exit(0);
    } else {
        std::cerr << "Pusher: " + fname + " successfully opened" << std::endl;
    }
    return myfile;
}
