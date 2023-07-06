/**
* \file    Tokamak.hpp
* \author  Xin Zhang, Princeton Plasma Physics Laboratory
* \date    December, 2020
*
* \brief   Implements the Tokamak class
* \details Stores the Tokamak plasma profile
*/
#include "Tokamak.hpp"

Tokamak::Tokamak(const std::string& dirIn, const std::string& eqdsk)
        : core_(nullptr), sol_(nullptr), mu_(2), Z_(1)
{
    core_ = new Core(dirIn, eqdsk);
    sol_ = new SOLPS(dirIn);
}

Tokamak::~Tokamak()
{
    if (core_) {
        delete core_;
    }
    if(sol_) {
        delete sol_;
    }
}

void Tokamak::setMainIon(int mu, int Z)
{
    mu_ = mu;
    Z_ = Z;
}

double Tokamak::getTemp(Vector pos, Particle background, size_t& init)
{
    double psiNow = (core_->machine_) -> getPsiNorm(pos);
    double rtn;
    if (psiNow <= INTERFACE){
        rtn = core_ -> getTemp(pos, background, init);
    } else {
        rtn = sol_ -> getTemp(pos, background, init);
    }
    return rtn;
}

double Tokamak::getDense(Vector pos, Particle background, size_t& init)
{
    double psiNow = (core_->machine_) -> getPsiNorm(pos);
    double rtn;
    if (psiNow <= INTERFACE){
        rtn = core_ -> getDense(pos, background, init);
    } else {
        rtn = sol_ -> getDense(pos, background, init);
    }
    return rtn;
}
