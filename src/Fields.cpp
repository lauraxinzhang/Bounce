/**
 * \file    Fields.cpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2020
 *
 * \brief   Implements the Fields class
 * 
 */

#include "Fields.hpp"

Vector Fields::getB(Vector &v)
{
    Vector rtn;
    return rtn;
}

bool Fields::isLimiter(Vector v)
{
    // TODO: implement this
    std::cerr << "isLimiter called in base class." << std::endl;
    return false;
}


void Fields::writeMat(std::string& filename, MatDoub& m)
{
    std::ofstream output;
    output.open(filename);
    
    
    if (!output.is_open()){
        std::cerr << "cannot open output file" << std::endl;
        exit(0);
    } else {
        output << m << std::endl;
    }
    output.close();
    return;
}

void Fields::writeVec(std::string& filename, VecDoub& v)
{
    std::ofstream output;
    output.open(filename);
    
    
    if (!output.is_open()){
        std::cerr << "cannot open output file" << std::endl;
        exit(0);
    } else {
        output << v << std::endl;
    }
    output.close();
    return;
}
