#
# Makefile for SOLFI
# g++ -o test -I/usr/pppl/gcc/9.3-pkgs/netcdf-cxx4-4.3.0/include test.cpp -L/usr/pppl/gcc/9.3-pkgs/netcdf-cxx4-4.3.0/lib -lnetcdf_c++4

# ----- Make Macros -----

CXX             = g++
CXXFLAGS        = -g -pedantic -w -Wall -Wextra -std=c++11 -O3

SRCDIR  = src/
BDDIR   = bin/

TARGETS = basic

EXECS   = $(BDDIR)basic.o

OBJECTS = $(BDDIR)Mesh.o $(BDDIR)TriagMesh.o $(BDDIR)QuadMesh.o $(BDDIR)UniMesh.o\
		  $(BDDIR)Particle.o $(BDDIR)Vector.o $(BDDIR)Matrix.o\
		  $(BDDIR)Containers.o $(BDDIR)Parser.o $(BDDIR)Fields.o $(BDDIR)FieldsTests.o\
		  $(BDDIR)GEquilibrium.o $(BDDIR)MeshTests.o\
		  $(BDDIR)Plasma.o $(BDDIR)SOL.o $(BDDIR)Core.o $(BDDIR)Tokamak.o $(BDDIR)Pusher.o

DEPS    =   $(SRCDIR)Mesh.hpp $(SRCDIR)TriagMesh.hpp $(SRCDIR)QuadMesh.hpp $(SRCDIR)UniMesh.hpp\
			$(SRCDIR)Particle.hpp $(SRCDIR)Parser.hpp\
			$(SRCDIR)Vector.hpp $(SRCDIR)Matrix.hpp $(SRCDIR)Containers.hpp \
			$(SRCDIR)Fields.hpp $(SRCDIR)FieldsTests.hpp $(SRCDIR)GEquilibrium.hpp\
			$(SRCDIR)MeshTests.hpp\
			$(SRCDIR)Plasma.hpp $(SRCDIR)SOL.hpp $(SRCDIR)Core.hpp\
			$(SRCDIR)Tokamak.hpp $(SRCDIR)Pusher.hpp\
			$(SRCDIR)delaunator.hpp $(SRCDIR)Constants.hpp $(SRCDIR)TestsUtil.hpp\

# ----- Make rules -----

all:	$(TARGETS)

clean:
	rm -rf $(TARGETS) $(OBJECTS) $(EXECS)

basic:	$(OBJECTS) $(BDDIR)basic.o
	$(CXX) $(CXXFLAGS) -o basic $(OBJECTS) $(BDDIR)basic.o


$(OBJECTS): $(BDDIR)%.o: $(SRCDIR)%.cpp $(DEPS)
	$(CXX) -c -o $@ $(CDFINC) $(CDFLIB) $< $(CXXFLAGS) $(LINKOPTS)

$(EXECS): $(BDDIR)%.o: $(SRCDIR)%.cpp $(DEPS)
	$(CXX) -c -o $@ $(CDFINC) $(CDFLIB) $< $(CXXFLAGS) $(LINKOPTS)

