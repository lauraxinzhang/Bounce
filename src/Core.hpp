/**
 * \file    Core.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    October, 2020
 *
 * \brief   Declares the Core class
 * \details Stores the core plasma profile
 */

#ifndef CORE_H_INCLUDED
#define CORE_H_INCLUDED 1

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
#include "Mesh.hpp"
#include "TriagMesh.hpp"
#include "QuadMesh.hpp"
#include "UniMesh.hpp"
#include "Fields.hpp"
#include "GEquilibrium.hpp"


class Core : public Plasma
{
public:
    
	Core(std::string dirIn  = "./input/", const std::string& eqdsk = "eqdsk");
    
    void setDirIn(std::string& dirIn);
    
    VecDoub readTrData(const char *  var);
    
    /**
     * \brief Get tempearture of the given species, at given flux cordinate
     * \param psiPol Normalized Poloidal Flux
     * \param background Sample particle of the background species
     */
    double getTemp(double psiPol, Particle background);
    
    /**
     * \brief Get density of the given species, at given flux cordinate
     * \param psiPol Normalized Poloidal Flux
     * \param background Sample particle of the background species
     */
    double getDense(double psiPol, Particle background);
    
    /**
     * \brief Overriding inherited  virtual function
     * \note t_init is not used in the Core class
     */
    double getTemp(Vector pos, Particle background, size_t& t_init);
    
    /**
    * \brief Overriding inherited  virtual function
     * \note t_init is not used in the Core class
     */
    double getDense(Vector pos, Particle background, size_t& t_init);
    
    void getValSep(double& Te, double& Ti, double& n);
    
    GEquilibrium * machine_;
    
    friend class Tokamak2P;
    
protected:
    std::string dirIn_;
    VecDoub trRMaj_;
    VecDoub trTe_;
    VecDoub trTi_;
    VecDoub trNe_;
    VecDoub trNi_;
    VecDoub trPsi_;
    
//    VecDoub trPsiPol_;
    XYMesh * mTe_;
    XYMesh * mTi_;
    XYMesh * mNe_;
    XYMesh * mNi_;
    
    
};

/**
 *\brief Housekeeping. Allows for exiting without memory leaks
 */
struct CoreException : public std::exception {
    int code_; ///< Exit code
    CoreException(int code): code_(code) {} // trivial constructor
    
    const char * what () const throw () {
        std::cerr << "CoreException encountered: " << std::endl;
        switch (code_){
            case 0:
                double zero(0);
        }
        return "Uncaught exceptions";
   }
};

#endif
