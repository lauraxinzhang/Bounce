/**
 * \file    SOL.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2020
 *
 * \brief   Declares the SOL class
 * \details Stores the SOL plasma profile
 */

#ifndef SOL_H_INCLUDED
#define SOL_H_INCLUDED 1

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

//class SOL
//{
//public:
//    SOL();
//
//	SOL(VecDoub& coords_, std::string dirOut  = "/Users/xinzhang/Documents/thesis/SOLFI/output/");
//
//    void setDirOut(std::string& dirOut);
//
//    void printMeshTriag();
//
//protected:
//    std::string dirOut_;
////    VecDoub coords_;
////    Mesh *mesh_;
////    VecDoub ne_, te_, nD_, ti_;
//};


class SOLPS : public Plasma
{
public:
    SOLPS(const std::string& dirIn);
    
    /**
     * \brief getting the values of the cell that the given mesh point is in
     * \return {ne, te, nD, ti}
     */
    VecDoub getVal(MeshPoint& p, size_t& t_init, bool histOut = false);
    
    /**
    * \brief Overriding inherited  virtual function
    */
    double getTemp(Vector pos, Particle background, size_t& t_init);
    
    /**
    * \brief Overriding inherited  virtual function
    */
    double getDense(Vector pos, Particle background, size_t& t_init);
    
    /**
     * \brief Print search history to file
     */
    void printHist(std::vector<size_t>& t_hist, std::vector<size_t>& q_hist);
    
    void printQuads();
    
    int getW();
    int getH();

private:
//    MatDoub ne_, te_, nD_, ti_;
    std::string dirOut_;
    int nx_, ny_;
    QuadMesh * quad_; // the inherited mesh_ is used instead
    TriagMesh * triag_;
    
    void setupTriag(std::vector<MatDoub>& crx, std::vector<MatDoub>& cry);
    
    
};


/**
 *\brief Housekeeping. Allows for exiting without memory leaks
 */
struct SOLException : public std::exception {
    int code_; ///< Exit code
    SOLException(int code): code_(code) {} // trivial constructor
    
    const char * what () const throw () {
        std::cerr << "SOLException encountered: " << std::endl;
        switch (code_){
            case 0: std::cerr << std::endl;
                exit(0); //normal exit
            case 1:
                std:: cerr << "Warning: Domain boundary encountered in search" << std::endl;
            case 2:
                std::cerr << "Fatal: point cannot be found within SOL domain." << std::endl;
                exit(2);
            case 3:
                std::cerr << "Warning: output files failed to open." << std::endl;
        
        }
        return "Uncaught exceptions";
   }
};
#endif
