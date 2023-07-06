/**
 * \file    Plasma.cpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    October, 2020
 *
 * \brief   Implements the Plasma class
 * \details Stores no private data; Provides methods for calculating
 *          collision frequencies for all children classes.
 */


#include "Plasma.hpp"

double Plasma::chandrasekharG(double x)
{
	double gauss = ((2 * x) / sqrt(PI)) * exp(-1 * pow(x, 2));
	return (erf(x) - gauss)/(2 * pow(x, 2));
}


void Plasma::xAndNu_ab(Particle& part, Particle& partB, double& xB, double& nu_ab)
{
	double dens, temp, ma, ea, mb, eb;
	ma = part.mass();
	ea = part.charge();
	mb = partB.mass();
	eb = partB.charge();

//	double rr = sqrt(part.pos().x() * part.pos().x() + part.pos().y() * part.pos().y() );
//	double zz = part.pos().z();
    
    
    dens = getDense(part.pos(), partB, part.searchInit_);
    temp = getTemp(part.pos(), partB, part.searchInit_) * QE; // convert to joules
    
	// First return value
	nu_ab =  dens * pow(ea, 2) * pow(eb, 2) * COULOMBLOG / (4 * PI * pow(EPSILON0, 2) * pow(ma, 2));

	// Second return value; temp is already in joules
	double vtB = sqrt(2 * temp / mb);
	xB = part.speed() / vtB;

	return;
}


double Plasma::nu_s(Particle& part, Particle& partB, double xB, double nu_ab)
{
	double ma, mb, temp;
//	double rr = sqrt(part.pos().x() * part.pos().x() + part.pos().y() * part.pos().y() );
//	double zz = part.pos().z();

    temp = getTemp(part.pos(), partB, part.searchInit_) * QE; // convert to joules

	ma = part.mass();
	mb = partB.mass();

	double g = chandrasekharG(xB);
	return g * (ma + mb) * nu_ab / (part.speed() * temp);
}

double Plasma::nu_D(Particle& part, double xB, double nu_ab)
{
	return nu_ab * (erf(xB) - chandrasekharG(xB)) / pow(part.speed(), 3);
}

double Plasma::nu_para(Particle& part, double xB, double nu_ab)
{
	return 2 * nu_ab * chandrasekharG(xB) / pow(part.speed(), 3);
}

double Plasma::root_D_para(Particle& part, double xB, double nu_ab)
{
    double ans = sqrt(2 * nu_ab * chandrasekharG(xB) / part.speed());
    return ans;
}

double Plasma::root_D_perp(Particle &part, double xB, double nu_ab)
{
    double ans = sqrt(nu_ab * (erf(xB) - chandrasekharG(xB)) / part.speed());
    return ans;
}

double Plasma::root_D_zero(double vtB, double nu_ab)
{
    double ans = sqrt(4 * nu_ab / (3 * sqrt(PI) * vtB));
    return ans;
}

void Plasma::thermalCoef(Particle part, double& nu_s_out, double& root_D_para_out, double& root_D_perp_out)
{
    Particle electron;
    Particle ion;
    ion.setSpec(2, 1);

    double xe, nu_e, xi, nu_i; // constants for electrons and main ions

    xAndNu_ab(part, electron, xe, nu_e);
    xAndNu_ab(part, ion, xi, nu_i);

    nu_s_out        = nu_s(part, electron, xe, nu_e) + \
                  nu_s(part, ion, xi, nu_i);
    root_D_para_out = root_D_para(part, xe, nu_e) + \
                  root_D_para(part, xi, nu_i);
    root_D_perp_out = root_D_perp(part, xe, nu_e) + \
                  root_D_perp(part, xi, nu_i);
    
//    if (psiNow <= INTERFACE){
//        core_ -> xAndNu_ab(part, electron, xe, nu_e);
//        core_ -> xAndNu_ab(part, ion, xi, nu_i);
//
//        nu_s        = core_ -> nu_s(part, electron, xe, nu_e) + \
//                      core_ -> nu_s(part, ion, xi, nu_i);
//        root_D_para = core_ -> root_D_para(part, xe, nu_e) + \
//                      core_ -> root_D_para(part, xi, nu_i);
//        root_D_perp = core_ -> root_D_perp(part, xe, nu_e) + \
//                      core_ -> root_D_perp(part, xi, nu_i);
//    } else {
//        sol_ -> xAndNu_ab(part, electron, xe, nu_e);
//        sol_ -> xAndNu_ab(part, ion, xi, nu_i);
//
//        nu_s        = sol_ -> nu_s(part, electron, xe, nu_e) + \
//                      sol_ -> nu_s(part, ion, xi, nu_i);
//        root_D_para = sol_ -> root_D_para(part, xe, nu_e) + \
//                      sol_ -> root_D_para(part, xi, nu_i);
//        root_D_perp = sol_ -> root_D_perp(part, xe, nu_e) + \
//                      sol_ -> root_D_perp(part, xi, nu_i);
//    }
    return;
}
