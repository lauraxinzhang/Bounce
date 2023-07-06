/**
 * \file    Plasma.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    October, 2020
 *
 * \brief   Declares the Plasma class
 * \details Stores no private data; Provides methods for calculating
 *          collision frequencies for all children classes.
 */

#ifndef PLASMA_H_INCLUDED
#define PLASMA_H_INCLUDED 1

#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "Particle.hpp"
#include "Vector.hpp"

class Plasma
{
public:
//    inline Plasma() : searchInit_(0){}
    
	/**
	 * \brief Template for getTemp.
	 * \note MUST be implemented by derived classes.
	 * \param pos 3D position vector of the query
	 * \param background An sample particle of the background species
	 * \param init (optional) cache for past searches
	 */
	inline virtual double getTemp(Vector pos, Particle background, size_t& init)
	{
		// throw exception and exit
		std::cerr << "unspecified plasma base class is referenced" << std::endl;
		return 0;
	}

	/**
	 * \brief Template for getDense.
	 * \note MUST be implemented by derived classes.
	 * \param pos 3D position vector of the query
	 * \param background An sample particle of the background species
	 * \param init (optional) cache for past searches
	 */
	inline virtual double getDense(Vector pos, Particle background, size_t& init)
	{
		// throw exception and exit
		std::cerr << "unspecified plasma base class is referenced" << std::endl;
		return 0;
	}

	/**
	 * \brief The chandrasekhar function
	 */
	double chandrasekharG(double x);


	/**
	 * \brief Calculates the constants xB and nu_ab
	 * \param part  Fast ion particle being pushed
	 * \param partB Background species
	 * \param xB   To be filled in; ratio of part.speed() and background vt
	 * \param nu_ab To be filled in; "base" collision frequency.
	 * \details Calculates xB and nu_ab, and update values of input references.
	 */
	void xAndNu_ab(Particle& part, Particle& partB, double& xB, double& nu_ab);

	/**
	 * \brief Calculates the current slowing down frequency of a given particle
	 * \param part  Particle that's currently being pushed
	 * \param partB Background particle; only need it for mass and charge
	 * \param xB    Ratio of particle speed and background thermal speed;
	 * \param nu_ab "base" collision frequency, basing on local density;
	 * \note        params xB and nu_ab should be calculated via Orbit::xAndNu_ab
	 */
	double nu_s(Particle& part, Particle& partB, double xB, double nu_ab);

	/**
	 * \brief Calculates perpendicular diffusion frequency
	 * \param part Particle currently being pushed
	 * \param xB    Ratio of particle speed and background thermal speed;
	 * \param nu_ab "base" collision frequency, basing on local density;
	 * \note        params xB and nu_ab should be calculated via Plasma::xAndNu_ab;;
     * Prefer to use root_D functions instead, to limit floating point errors.
	 */
	double nu_D(Particle& part, double xB, double nu_ab);

	/**
	 * \brief Calculates parallel diffusion frequency
	 * \param part Particle currently being pushed
	 * \param xB    Ratio of particle speed and background thermal speed;
	 * \param nu_ab "base" collision frequency, basing on local density;
	 * \note        params xB and nu_ab should be calculated via Plasma::xAndNu_ab
	 */
	double nu_para(Particle& part, double xB, double nu_ab);
    
    /**
     * \brief Square-root of the perpendicular diffusion coefficient
     */
    double root_D_perp(Particle& part, double xB, double nu_ab);
    
    /**
    * \brief Square-root of the parallel diffusion coefficient
    */
    double root_D_para(Particle& part, double xB, double nu_ab);
    
    /**
     * \brief Square-root of the diffusion coefficient at v/vt << 1 limit
     * \note Ensures that the diffusion coefficients are well behaved as v -> 0
     */
    double root_D_zero(double vtB, double nu_ab);
    
    /**
     *\brief Calculates the coefficients for thermal collisions
     *\details Calculates the combined thermal coefficients from electrons and ions. Ions are assumed to be D+ in the base class implementation. Inherited classes should override this to change the default behavior.
     * \param part Test particle - the velcoty and position are needed.
     * \param nu_s_out Slowing down frequency. Reference to pre-initialized double. Value modified after function call.
     * \param root_D_para_out Square root of D_parallel. Modified after function call.
     * \param root_D_perp_out Square root of D_perp. Modified after function call.
     */
    virtual void thermalCoef(Particle part, double& nu_s_out, double& root_D_para_out, double& root_D_perp_out);
};

#endif
