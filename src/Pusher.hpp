/**
 * \file    Pusher.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    October, 2020
 *
 * \brief   Declares the Pusher class
 */

#ifndef Pusher_H_INCLUDED
#define Pusher_H_INCLUDED 1

#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "Plasma.hpp"
#include "Core.hpp"
#include "SOL.hpp"
#include "Fields.hpp"
#include "Particle.hpp"
#include "Parser.hpp"
#include "Tokamak.hpp"


class Pusher
{
public:
    /**
     *\brief Minimal contructor for particle pushing. Plasma_ is set to null: no thermal plasma provided.
     *\param fields Magnetic / electric field geometry
     */
    Pusher(Fields& fields);

    /**
     * \brief Constrcutor for a general configuration of Plasma and Fields
     * \param plasma The thermal plasma. Run-time polymorphic - see documentation for the Plasma class
     * \param fields Magnetic / electric field geometry
     */
    Pusher(Plasma& plasma, Fields& fields);
    
    /**
     * \brief Constructor for the core-sol coupled plasma background. Short hand for me
     * \param dirIn Input directory to look for plasma profiles for sol and core
     * \param eqdsk Name of the equilibrium file, in input directory
     */
	Pusher(const std::string& dirIn, const std::string& eqdsk);
    
//    ~Pusher();
    
    /**
     * \brief Change output directory
     */
    void setDirOut(std::string out);
    
    
    /**
     * \brief Initializes the particle list with Gaussian velocity
     * \param nparts Number of particles
     * \param energy Energy of the particles (eV)
     * \param pos         Initial position of particles (shared by all particles)
     * \param mu           Mass number, 0 for electrons
     * \param Z             Charge number, -1 for electrons
     */
    void initParts(int nparts, double energy, Vector pos, int mu, int Z);
    
    /**
     * \brief Pushes a single particle, returns final position
     * \param part The particle to be pushed
     * \param iter Number of iterations
     * \param collide Whether the particle interacts with thermal background
     * \param write Whether the trajectory should be written to file
     */
    void pushSingle(Particle& part, int iter,  bool collide, std::string prefix, bool write = 1);
    
    /**
     *\brief Pushes all particles in the partList_;
     *TODO: not implemented
     */
    void pushAll();
    
//    /**
//     *\brief Calculates the coefficients for thermal collisions
//     *\note specific to the core + sol plasma structure. Uses the INTERFACE value defined in constants.h
//     *\details Looks for plasma background in core_ if the psi_pol of the particle position is less than INTERFACE; looks
//     * for sol_ plasma profile otherwise. Contributions from electrons and main ion are added.
//     * \param part Test particle - the velcoty and position are needed.
//     * \param nu_s Slowing down frequency. Reference to pre-initialized double. Value modified after function call.
//     * \param root_D_para Square root of D_parallel. Modified after function call.
//     * \param root_D_perp Square root of D_perp. Modified after function call.
//     * TODO: pull the core and edge switch to a different file. Keep this class general for any plasma configuration.
//     */
//    void thermalCoef(Particle part, double& nu_s, double& root_D_para, double& root_D_perp);
    
    /**
     *\brief Safely open a new file for outputting
     *\param fname Name of the output file. File is located in output directory dirOut_;
     *\note  Throws an error and exit the program if file cannot be opened.
     */
    std::ofstream safeOpen(std::string fname);
    

private:
	Fields * fields_;
    Plasma * plasma_;

    std::vector<Particle> partList_;
    
    std::string dirOut_;
    
};


#endif
