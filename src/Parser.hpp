/**
 * \file    Parser.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2020
 *
 * \brief   Declares the Parser class
 * \details Handles file input parsing. 
 *			Transfers file contents to data containers only, all subsequent
 *          processing happen in other classes.
 */

#ifndef PARSER_H_INCLUDED
#define PARSER_H_INCLUDED 1

#include "Containers.hpp"
#include <cstdio>
#include <time.h> // for seeding the random engine.
#include <random>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

class Parser
{
public:

    /**
     *\brief Stores information about domain boundary and limiter coordinates
     */
    struct limiter
    {
        int nbbbs, limitr;
        VecDoub rbbbs_, zbbbs_, rlim_, zlim_;
        
        /**
         * \brief Constructor for limiter structure
         * \note All vectors initialized to desired size, but not filled in until
         *       Parser::parseEqdsk is called.
         */
        inline limiter(int nbbbs, int limitr)
                    :nbbbs(nbbbs), limitr(limitr),
                    rbbbs_(nbbbs), zbbbs_(nbbbs),
                    rlim_(limitr), zlim_(limitr) {}
        
        /**
         *\brief Default constructor for limiter structure.
         *\note Needed for default initialization of the eqdsk struct
         */
        inline limiter() :nbbbs(0), limitr(0){}
    };
    
	/**
	 * \brief Stores all informations from an eqdsk file
	 */
	struct eqdsk
	{
		int nw, nh;
		double rdim, zdim, rcentr, rleft, zmid;
		double rmaxis, zmaxis, simag, sibry, bcentr;
		double current, simag1, xdum, rmaxis1, xdum1;
		double zmaxis1, xdum2, sibry1, xdum3, xdum4;

        // all variables with trailing underscores are vectors
		VecDoub fpol_, siGrid_, pres_, ffprim_, pprime_;
		VecDoub rGrid_, zGrid_, qpsi_;
		MatDoub psirz_;

        /**
         * \brief data structure that constains boundary and limiter coordinates
         */
        limiter limit_;
        
//        /**
//         *\brief Default constructor
//         */
//        inline eqdsk()
//                :nw(0), nh(0),
//                rdim(0), zdim(0), rcentr(0), rleft(0), zmid(0),
//                rmaxis(0), zmaxis(0), simag(0), sibry(0), bcentr(0),
//                current(0), simag1(0), xdum(0), rmaxis1(0), xdum1(0),
//                zmaxis1(0), xdum2(0), sibry1(0), xdum3(0), xdum4(0),
//                fpol_(0), siGrid_(0), pres_(0), ffprim_(0), pprime_(0),
//                rGrid_(0), zGrid_(0), qpsi_(0), psirz_(0, 0) {}

        /**
         * \brief Constructor for eqdsk structure
         * \note All vectors initialized to desired size, but not filled in until
         *       Parser::parseEqdsk is called.
         */
		inline eqdsk(int nw, int nh)
                :nw(nw), nh(nh),
                rdim(0), zdim(0), rcentr(0), rleft(0), zmid(0),
                rmaxis(0), zmaxis(0), simag(0), sibry(0), bcentr(0),
                current(0), simag1(0), xdum(0), rmaxis1(0), xdum1(0),
                zmaxis1(0), xdum2(0), sibry1(0), xdum3(0), xdum4(0),
                fpol_(nw), siGrid_(nw), pres_(nw), ffprim_(nw), pprime_(nw),
                rGrid_(nw), zGrid_(nh), qpsi_(nw), psirz_(nw, nh)
        {
            return;
        }
        
        
        /**
         * \brief Access for the limiter structure parsed from eqdsk
         * \note Same as directly calling self.limit_
         */
        inline limiter& getLimiter()
        {
            return limit_;
        }
        
    };
    

    /**
     * \brief Default constuctor
     * \note Does nothing, since Parser class doesn't have private members.
     */
	Parser();
    
    std::ifstream safeOpen(const std::string& dirIn, const std::string& fname);
    
    std::ifstream safeOpen(const std::string& path);

    /**
     * \brief Parses an input eqdsk file, returns an eqdsk structure
     */
	eqdsk parseEqdsk(const std::string& path);

    /**
     * \brief Parses an .dat file produced from SOLPS
     * \note  All matrices are flattened
     */
    VecDoub parseDat(const std::string& dirIn, const std::string& fname);
    
    /**
     *\brief Parses an .dat file produced from SOLPS
     *\note Values are stored in a MatDoub. unflattened.
     */
    MatDoub parseDatMat(const std::string& dirIn, const std::string& fname);
    
    VecDoub parseVec(const std::string& dirIn, const std::string& fname);

private:
//	std::string magPath_;
//	std::string limitPath_;

};

#endif
