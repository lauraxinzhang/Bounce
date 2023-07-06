/**
 * \file    Containers.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    June, 2020
 *
 * \brief   Implements the VecDoub and MatDoub classes
 * 
 */
#include "Containers.hpp"

MatDoub::MatDoub(int nw, int nh, double a)
		:nw_(nw), nh_(nh)
{
    // Row major
	data_.resize(nh_);
	std::vector<double> rows(nw_, a);
	for (int i = 0; i < nh_; i++) {
		data_.at(i).resize(nw_);
		data_.at(i) = rows;
	}
	return;
}

MatDoub::MatDoub(VecDoub& vec)
        :nw_(vec.size()), nh_(1)
{
    data_.resize(1);
    data_.at(0).resize(nw_);
    data_.at(0) = vec;
    return;
}

double& MatDoub::at(int w, int h)
{
	return data_.at(h).at(w);
}

VecDoub& MatDoub::getRow(int h)
{
	return data_.at(h);
}

VecDoub MatDoub::getCol(int w)
{
    VecDoub rtn(nh_);
    for (int i = 0; i < nh_; i++){
        rtn.at(i) = at(w, i);
    }
    return rtn;
}

int MatDoub::getW()
{
	return nw_;
}

int MatDoub::getH()
{
	return nh_;
}

void MatDoub::resize(int nw, int nh)
{
    data_.resize(nh);
    for (int i = 0; i < nh; i++) {
        data_.at(i).resize(nw);
    }
    nw_ = nw;
    nh_ = nh;
}

VecDoub MatDoub::dot(VecDoub& right)
{
    int size = right.size();
    try {
        if (size != nw_){
            throw CtnrException(0);
        }
    } catch (CtnrException& e) {
        e.what();
    }
    VecDoub rtn(nh_);
    for (int i=0; i < nh_; i++){
        double sum(0);
        VecDoub row = getRow(i);
        for (int j = 0; j< nw_; j++){
            sum += row.at(j) * right.at(j);
        }
        rtn.at(i) = sum;
    }
    return rtn;
}
