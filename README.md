# Bounce [![DOI](https://zenodo.org/badge/278388974.svg)](https://zenodo.org/badge/latestdoi/278388974)

Bounce bounce bounce a ball, all the bananas!


Particle full orbit tracer in full machine 3D geometry.

## Compiling Bounce

Requires a g++ compiler that supports the c++11 standard. No external libraries are needed. Run the following command to make all targets:

```bash
make
```
This will compile bother targets ```basic``` and ```buildsol```. 

To clean up build objects, do
```bash
make clean
```

## Usage

### basic
```basic``` executable constructs a particle pusher object that includes a full tokamak model. The particle is then evolved in time until maximum iteration is hit or if the particle is lost to the limiter.

Update particle initializations in basic.cpp. Recompile the program and run:

```bash
make
./basic
```

### API
To use the provided Plasma and Fields user classes, first modify the main controller file, in basic.cpp. 

To implement custom configurations, inherit Plasma and Fields class (run-time polymorphic), and override the virtual functions. Everything else should work as usual.


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
=======
