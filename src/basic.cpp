#include <cstdio>
#include <time.h> // for seeding the random engine.
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "Containers.hpp"
#include "Particle.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "Parser.hpp"
#include "Fields.hpp"
#include "FieldsTests.hpp"
#include "Pusher.hpp"
#include "Plasma.hpp"


void runTests()
{
    FieldsTests ftest;
    ftest.runAll();
}


int main() {
    runTests();

    std::string dirIn("./input/");
    std::string dirOut("./output/");
  
    std::string nstu = "gD3D.geqdsk"; // ok

    // Setup the particle pusher driver
    GEquilibrium gfile(dirIn + nstu); // construct magnetic field
    Pusher pusher(gfile);
    pusher.setDirOut(dirOut);

    // Initialize a particle
    Vector pos(1.5, 0, 0); // change this

    Particle part(pos, pos, 1, 1); // position, velocity, mass, charge

    double energy = 4000;
    double vt = sqrt(2 * energy * QE / part.mass());
    Vector vel(vt, 0, 0);
    part.setVel(vel);

    // push the particle
    int maxiter = 100000;
    bool collision = false; // lemme know if at some point you want collisional 
    pusher.pushSingle(part, maxiter, collision, "example", 1);
}


