/**
 * \file    Tokamak.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    December, 2020
 *
 * \brief   Declares the Tokamak class
 * \details Stores the Tokamak plasma profile
 */

#ifndef TOKAMAK_H_INCLUDED
#define TOKAMAK_H_INCLUDED 1

#include <cstdio>
#include <time.h> // for seeding the random engine.
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <assert.h>

#include "Plasma.hpp"
#include "Containers.hpp"
#include "Parser.hpp"
#include "Constants.hpp"
#include "SOL.hpp"
#include "Core.hpp"

class Tokamak : public Plasma
{
public:
    Tokamak(const std::string& dirIn, const std::string& eqdsk);
    
    ~Tokamak();
    
    void setMainIon(int mu, int Z);
    
    double getTemp(Vector pos, Particle background, size_t& init);
    
    double getDense(Vector pos, Particle background, size_t& init);
    
//    void thermalCoef(Particle part, double& nu_s_out, double& root_D_para_out, double& root_D_perp_out);
    
private:
    Core * core_;
    SOLPS * sol_;
    
    int mu_; // mass number of main ion
    int Z_; // charge number of main ion
};




#endif
