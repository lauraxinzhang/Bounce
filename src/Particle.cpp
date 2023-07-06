/**
 * \file    Particle.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2017
 *
 * \brief   Implements the Particle class
 * 
 */

#include "Particle.hpp"

Particle::Particle()
		:pos_(), vel_(), mu_(0), dt_(-1),
        lost_(false), searchInit_(0), generator_(int(time(NULL)))
{
	charge_ = -1 * QE;
	mass_ = ME;
} 

Particle::Particle(const Vector& pos, const Vector& vel, int mu, int Z)
		:charge_(Z * QE),
		 pos_(pos.x(), pos.y(), pos.z()),
		 vel_(vel.x(), vel.y(), vel.z()),
		 lost_(false), mu_(mu), dt_(-1),
         searchInit_(0), generator_(int(time(NULL)))
{
    if (mu == 0){
        mass_ = ME;
    } else {
        mass_ = mu * MI;
    }
}

Vector Particle::pos() const
{
	return pos_;
}

Vector Particle::vel() const
{
	return vel_;
}

//bool Particle::spec() const
//{
//	return species_;
//}

double Particle::mass() const
{
	return mass_;
}

double Particle::charge() const
{
	return charge_;
}

double Particle::energy() const
{
	double result = 0.5 * mass() * pow(speed() , 2)/QE;
	return result;
}

double Particle::setDt(double Btypical)
{
    double fLamor;
    if (mu_ == 0){ // if particle is an electron
        fLamor = 2.8E6 * (Btypical * TESLATOGAUSS);
    } else {
        fLamor = 1520 / mu_ * (Btypical * TESLATOGAUSS);
    }
    // constants from NRL p28 ( Btypical is converted to Gauss.
    double TLamor = 1.0/fLamor;
    double dt = TLamor / NPERORBIT;
    dt_ = dt;
    return dt;
}

double Particle::speed() const
{
	return vel().mod();
}

void Particle::setPos(const Vector& right)
{
	pos_.setX(right.x());
	pos_.setY(right.y());
	pos_.setZ(right.z());
	return;
}

void Particle::setVel(const Vector& right)
{
	vel_.setX(right.x());
	vel_.setY(right.y());
	vel_.setZ(right.z());
	return;
}

void Particle::setSpec(int mu, int Z)
{
    if (mu != 0){
        mu_ = mu;
        mass_ = mu * MI;
    } else {
        mass_ = ME;
    }
    charge_ = Z * QE;
	return;
}

void Particle::lost()
{
	lost_ = true;
	return;
}

bool Particle::isLost()
{
	return lost_;
}

double Particle::magMoment(Vector& BField) const
{
	Vector vperp = vel().perp(BField);
	double perp = vperp.mod();
	double B = BField.mod();

	double mu = (mass() * perp * perp)/ (2 * B);

	return mu;
}

bool Particle::operator==(const Particle& right)
{
	bool posB = (pos() == right.pos());
	bool velB = (vel() == right.vel());
	bool specB = (mass() == right.mass() && charge() == right.charge());
	return (posB && velB && specB);
}

void Particle::move(Vector& E, Vector& B)
{
    if (fabs(dt_ + 1) < PRECISION){ // dt_ has not been initialized
        double B0 = B.mod();
        setDt(B0);
    }
	double qprime;
	
	qprime = dt_ * charge_ / (2.0 * mass_);

	Vector h = B * qprime;
	// std::cerr << "h:" << h << std::endl;

	double hMod = h.mod();
	// std::cerr << "hMod:" << hMod << std::endl;

	Vector s = (h * 2.0)/(1.0 + pow(hMod,2));
	// std::cerr << "s:" << s << std::endl;

	Vector u = vel_ + E * qprime;
	Vector uuh = u + (u.cross(h));
	// std::cerr << "uuh:" << uuh << std::endl;

	Vector uPrime = u + uuh.cross(s);
	// std::cerr << "u':" << uPrime << std::endl;

	vel_ = uPrime + E * qprime;
	pos_ = pos_ + vel_ * dt_; // used updated velocity
	// std::cerr << pos_ << std::endl;
	return;
}

/* THIS IS WRONG - DIFFERENTIAL GEOMETRY WAS NOT HANDLED CORRECTLY

void Particle::moveCyl(Vector& E, Vector& B, const double dt)
{

	double oldR     = pos_.x();
	double oldTheta = pos_.y();
	double oldZ     = pos_.z();

	Vector newpos(oldR, 0, oldZ);

	setPos(newpos);
	move(E, B, dt);

	double newX = pos_.x();
	double newY = pos_.y();
	double newZ = pos_.z();

	// convert pos vector back to cylindrical coordinates

	double newR   = sqrt(pow(newX, 2) + pow(newY, 2));
	// changed from oldR to newX in the following line - 11/14/17

	double dTheta = atan( newY / newX );

	Vector newPos(newR, oldTheta + dTheta, newZ);
	// std::cerr << newPos << std::endl;

	// rotate velocity vector to align with current reference frame:
	double ux = vel_.x();
	double uy = vel_.y();
	double uz = vel_.z();

	double vx = cos(dTheta) * ux - sin(dTheta) * uy;
	double vy = sin(dTheta) * ux + cos(dTheta) * uy;


	Vector newVel(vx, vy, uz);

	setVel(newVel);
	setPos(newPos);

	return;
}

*/

void Particle::moveCyl(Vector& E, Vector& B)
{
	/* The particle pusher is now entirely in cartesian coordinates. No
	artificial rotations of any kind. Toroidal motions should already be
	resolved in the 3D Boris algorithm. All fields (accelerations) are 
	converted to Cartesian before applied to the pusher.
	*/

	// Now assume that particle position and velocity vectors are in xyz

	// Magnetic field is still in R, phi, Z
	// Electric field is solved in RZ plane, forcing rotational
	// symmetry, therefore no field in phi direction. We'll keep electic
	// field in R, phi, Z coordinates too.

   	Vector Enew;
   	Vector Bnew;

   	E.cyl2Cart(pos(), Enew);
   	B.cyl2Cart(pos(), Bnew);
   	
   	// The good old Boris push.
   	move(Enew, Bnew);

   	return;
}

void Particle::scatter(Vector& B)
{
    std::cerr << "Warning: Large angle collision is deprecated. It is NOT gauranteed to work correctly" << std::endl;
	std::bernoulli_distribution distribution(0.5);
	bool sign = distribution(generator_); // choose the orientation of the new vector by random;
	
	Vector axis = vel().cross(B).normalize();

	Vector newv = vel().turn(axis, sign);
	setVel(newv);
	return;
}

void Particle::collide(double nu_s, double root_D_para, double root_D_perp)
{
    assert(dt_ != -1); // make sure that dt has been initialized
    // update particle velocity: slowing down
    Vector vslow = vel_ * nu_s * dt_;
    
    Matrix identity;
    identity.diagonal(1, 1, 1);

    double v = speed();
//    std::cerr << "speed: " << speed() << "vslow: " << vslow.mod() << std::endl;
    Matrix vv = vel_.tensor(vel_);
    Matrix parallel =  vv * pow(v, -2);
    Matrix perp     = identity - parallel;
    Matrix diffusion = parallel * root_D_para + perp * root_D_perp;

    std::normal_distribution<double> distribution(0, sqrt(dt_)); // generate a Gaussian distributed velocity
    double wx = distribution(generator_); // generate 3 normal distributed velocities.
    double wy =  distribution(generator_);
    double wz = distribution(generator_);
    Vector wiener(wx, wy, wz);

    Vector diffused = diffusion.dot(wiener);

    Vector updated = vel_ - vslow + diffused;
    setVel(updated);
    return;
}
