/**
 * \file    Parser.cpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2020
 *
 * \brief   Implements the Parser class
 * 
 */
#include "Parser.hpp"

Parser::Parser()
{
    return;
}

Parser::eqdsk Parser::parseEqdsk(const std::string& path)
{
    std::ifstream myfile = safeOpen(path);

    std::string a;
    // dimension information
    int nw, nh;
    
    getline(myfile, a);
    std::stringstream topRow(a);
    std::string token;
    std::vector<std::string> tokens;
    while(getline(topRow, token, ' '))
    {
        if (!token.empty()){
            tokens.push_back(token);
        }
    }
    size_t len = tokens.size();
    try {
        nw = std::stoi(tokens.at(len - 2));
        nh = std::stoi(tokens.at(len - 1));
    } catch (std::exception& e) {
        std::cerr << "Exception: " << e.what();
        std::cerr << "nw and nh information cannot be found in the eqdsk." << std::endl;
        std::cerr << "Please check equilibrium file format and try again." << std::endl;
        exit(0);
    }
    
//    std::cout << "New parser: nw, nh " << nw << ',' << nh << std::endl;

    // create Parser object
	Parser::eqdsk rtn(nw, nh);
    
    myfile >> std::setprecision(10) >> rtn.rdim >> rtn.zdim;
    myfile >> rtn.rcentr >> rtn.rleft >> rtn.zmid;

    myfile >> rtn.rmaxis >> rtn.zmaxis >> rtn.simag >> rtn.sibry >> rtn.bcentr;
    
    myfile >> rtn.current >> rtn.simag1 >> rtn.xdum >> rtn.rmaxis1 >> rtn.xdum1;
    myfile >> rtn.zmaxis1 >> rtn.xdum2 >> rtn.sibry1 >> rtn.xdum3 >> rtn.xdum4;
    
//    std::cout << std::scientific << "rleft: " << rtn.rleft << " zmid:" << rtn.zmid << std::endl;
//    std::cout << "rmaxis: " << rtn.rmaxis << " zmaxis:" << rtn.zmaxis << std::endl;
    
    double dsi = (rtn.sibry - rtn.simag) / (rtn.nw - 1);

    for (int i = 0; i < nw; ++i){
        double d;
        myfile >> std::setprecision(10) >> d;
        rtn.fpol_.at(i) = d;
        rtn.siGrid_.at(i) = rtn.simag + i * dsi;
    }

    for (int i = 0; i < nw; ++i){
       double d;
       myfile >> std::setprecision(10) >> d;
       rtn.pres_.at(i) = d;
    }

    for (int i = 0; i < nw; ++i){
       double d;
       myfile >> std::setprecision(10) >> d;
       rtn.ffprim_.at(i) = d;
    }

    for (int i = 0; i < nw; ++i){
       double d;
       myfile >> std::setprecision(10) >> d;
       rtn.pprime_.at(i) = d;
    }

    for (int j = 0; j < nh; ++j){
        for (int i = 0; i < nw; ++i){
            double d;
            myfile >> std::setprecision(10) >> d;
            rtn.psirz_.at(i, j) = d;
        }
    }

    double dr = rtn.rdim / (rtn.nw - 1);
    for (int i = 0; i < nw; ++i){     
        rtn.rGrid_.at(i) = rtn.rleft + i * dr;
    }

    double dz = rtn.zdim / (rtn.nh - 1);
    for (int j = 0; j < nh; ++j){ 
       rtn.zGrid_.at(j) = (rtn.zmid - rtn.zdim/2.0) + j * dz;
    }
    
    VecDoub qpsi(nw);
    for (int i = 0; i < nw; ++i){
        double d;
        myfile >> std::setprecision(10) >> d;
        rtn.qpsi_.at(i) = d;
    }

    // number of boundary points and limiter points'
    int nbbbs, limitr;
    myfile >> std::setprecision(10) >> nbbbs >> limitr;
    limiter limit(nbbbs, limitr);
    
    for (int i = 0; i < limit.nbbbs; ++i){
        double d;
        myfile >> std::setprecision(10) >> d;
        limit.rbbbs_.at(i) = d;
        
        myfile >> std::setprecision(10) >> d;
        limit.zbbbs_.at(i) = d;
    }

    for (int i = 0; i < limit.limitr; ++i){
        double d;
        myfile >> std::setprecision(10) >> d;
        limit.rlim_.at(i) = d;
        myfile >> std::setprecision(10) >> d;
        limit.zlim_.at(i) = d;
    }

    rtn.limit_ = limit;
    
    std::cerr << "eqdsk file: " << path << " successfully parsed!" << std::endl;
	return rtn;
}

VecDoub Parser::parseDat(const std::string& dirIn, const std::string& fname)
{
    std::ifstream myfile = safeOpen(dirIn, fname);

    // grab the first line, containing indices
    std::string line;
    std::getline( myfile, line );
    std::istringstream is( line );

    int iw(0), nw(0);
    while (is >> iw){
        nw = iw; // keep rewriting the value of nw until done
    }

    VecDoub val(nw);
    int nh(-5); // an invalid starting value
    // Parse through the rest of the file
    while ( std::getline( myfile, line ) ){
        std::istringstream is( line );
        int ih; // index of the row
        is >> ih;
        if (nh == -5){ // nh is uninitialized
            nh = ih;
            val.resize(nh * nw);
        }
        if (ih != nh && ih != -1){ // if we're not at the first or last row
            double d;
            is >> d;
            for (int iw = 2; iw < nw + 2; ++iw){
                is >> d;
                int index = ih * nw + iw - 2;
                val.at(index) = d;
            }
        } // otherwise do nothing
    }
    return val;
}

MatDoub Parser::parseDatMat(const std::string& dirIn, const std::string& fname)
{
    std::ifstream myfile = safeOpen(dirIn, fname);

    // grab the first line, containing indices
    std::string line;
    std::getline( myfile, line );
    std::istringstream is( line );

    int iw(0), nw(0);
    while (is >> iw){
        nw = iw; // keep rewriting the value of nw until done
    }
    
    MatDoub rtn(nw, 0);
    int nh(-5); // an invalid starting value
    // Parse through the rest of the file
    while ( std::getline( myfile, line )){
        if (! line.empty()){
            std::istringstream is( line );
            
            int ih; // index of the row
            is >> ih;
            if (nh == -5){ // nh is uninitialized
                nh = ih;
                rtn.resize(nw, nh);
            }
            if (ih != nh && ih != -1){ // if we're not at the first or last row
                double d;
                is >> d;
                for (int iw = 2; iw < nw + 2; ++iw){
                    is >> d;
                    rtn.at(iw - 2, ih) = d;
                }
            } // otherwise, row is guard cell; do nothing.
        }
    }
    return rtn;
}

VecDoub Parser::parseVec(const std::string& dirIn, const std::string& fname)
{
    std::ifstream myfile = safeOpen(dirIn, fname);
    VecDoub rtn(0);
    double d;
    while (myfile >> d){
        rtn.push_back(d);
    }
    return rtn;
}

std::ifstream Parser::safeOpen(const std::string& dirIn, const std::string& fname)
{
    std::ifstream myfile;
    std::string path = dirIn + fname;
    myfile.open(path); // open eqdsk file

    if (!myfile.is_open()) {
        std::cerr << "Unable to open file: " + fname << std::endl;
        exit(0);
    } else {
//        std::cerr << "Parser: " + fname + " successfully loaded" << std::endl;
    }
    return myfile;
}

std::ifstream Parser::safeOpen(const std::string& path)
{
    std::ifstream myfile;
    myfile.open(path); // open eqdsk file

    if (!myfile.is_open()) {
        std::cerr << "Unable to open file: " + path << std::endl;
        exit(0);
    } else {
//        std::cerr << "Parser: " + path + " successfully loaded" << std::endl;
    }
    return myfile;
}
