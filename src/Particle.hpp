/**
 * \file    Particle.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2017
 *
 * \brief   Declares the Particle class
 * 
 */

#ifndef PARTICLE_H_INCLUDED
#define PARTICLE_H_INCLUDED 1

#include <iostream>
#include <time.h> // for seeding the random engine.
#include <random>
#include "Vector.hpp"
#include "Matrix.hpp"
#include "Constants.hpp"
#include <assert.h>

/**
 * \brief  A class for individual particle objects
 */
class Particle
{
    public:

    	/**
    	 * \brief Default constructor. 
    	 * \details Calls the default constructor for Vector class to initialize both
    	 *          pos and vel to [0,0,0]. Species default to electron.
    	 */
    	Particle();   // default constructor
    	// ~Particle();

        /**
         * \brief Another constructor for particles other than proton and electrons
         * \param pos    position of particle
         * \param vel    velocity of particle
         * \param mass   Mass number of particle (in proton masses)
         * \param charge Charge number of particle
         *
         */
        Particle(const Vector& pos, const Vector& vel, int mass, int charge);

    	/**
    	 * \brief Data getter for position
    	 */
    	Vector pos() const;

    	/**
    	 * \brief Data getter for velocity
    	 */
    	Vector vel() const;

        /**
         * \brief Data getter for mass
         */
        double mass() const;

        /**
         * \brief Data getter for charge
         */
    	double charge() const;

        /**
         * \brief Calculates particle energy (in eV)
         */
    	double energy() const;

        /**
         * \brief Calculates time step for moving
         * \param Btypical a "typical" value of the magnetic field to calculate Lamor frequency with.
         */
        double setDt(double Btypical);
    
        /**
         * \brief Calculated particle speed
         */
        double speed() const;

    	void setPos(const Vector& right);
    	void setVel(const Vector& right);
    
        /**
         *\brief Sets the species of the particle
         *\param mu Mass number of the particle. 0 for electrons, otherwise the atomic mass number
         *\param Z Charge of the particle. -1 for electrons
         */
    	void setSpec(int mu, int Z);

        /**
         * \brief Set the lost_ flag to true if the particle is lost.
         */
        void lost();

        /**
         * \brief Return whether the particle is lost
         */
        bool isLost();

        /**
         * \brief Calculate 0th order magnetic moment at the given B field
         * \note  B Field value must be the field at the current position
         */
        double magMoment(Vector& BField) const;

        /**
         * \brief Comparison operator. Compares pos, vel, and spec
         * \note  state of random number generator is not compared
         */
    	bool operator==(const Particle& right);

        /**
         * \brief Evolve particle position and velocity
         * \param E local electric field (SI units)
         * \param B local magnetic field (SI units)
         * \note  Everything is handled in cartesian geometry. If E and B
         *        are defined in toroidal geometry, the field vectors need
         *        to be converted to cartesian before this method can be called.
         */
    	void move(Vector& E, Vector& B);

        /**
         * \brief Particle pusher in cylindrical geometry
         * \param E local electric field (Toroidal, SI)
         * \param B local magnetic field (Toroidal, SI)
         * \note  Cast (R, phi, Z) field vectors into cartesian geometry according
         *        to the current location of the particle, then call Particle::move
         */
        void moveCyl(Vector& E, Vector& B);

        /**
         * \brief Turn the particle pitch angle by 90 degrees
         * \param B Local magnetic field
         * \note  DEPRECATED. Coulomb collision is now handled directly with
         *        Wiener processes.
         */
        void scatter(Vector& B);
    
        /**
         * \brief a small angle scattering, described by the given 3 coefficients
         * \param nu_s Slowing down frequency
         * \param root_D_para Parallel diffusion coefficient
         * \param root_D_perp Perpendicular diffusion coefficient
         */
        void collide(double nu_s, double root_D_para, double root_D_perp);

        /**
        * \brief Cache for search (used in unstructured grids)
        */
        size_t searchInit_;
    
    private:
        double   mass_;
        int      mu_; // mass number
        double   charge_;

    	Vector pos_; // pos at time step k
    	Vector vel_; // vel at time step k-1/2
//    	bool   species_;
    
        double dt_;

        bool   lost_; // always start with false. Change to true when limiter is reached.
    
        std::default_random_engine generator_; // a random number generator, seeded at Particle construction;
};

#endif // PARTICLE_H_INCLUDED
