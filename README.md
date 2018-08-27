# DGSWEM V2
[![CircleCI](https://circleci.com/gh/UT-CHG/dgswemv2.svg?style=svg&circle-token=0479b7746d69a87e977936dd4b6105be5b2e6316)](https://circleci.com/gh/UT-CHG/dgswemv2)

## Getting Started

dgswem-v2 has a considerable number of dependencies. If you are simply interested in giving dgswemv2 a quick test. We recommend one of two options. First, if one does not want to build all of the dependencies outlined below, running
```sh
cmake -DCMAKE_PREFIX_PATH=<YAML_CPP_INSTALL_PATH>\
      -DBUILD_EXAMPLES=On\
    <DGSWEMV2_ROOT>
```
will build the code to run using a Runge-Kutta discontinuous Galerkin discretization in serial along with the manufactured solution.

Another alternative is to use our Docker image, which contains the necessary dependencies to build all dgswemv2 targets. Simply, run
```sh
docker pull bremerm31/dgswemv2:latest
docker run -it bremerm31/dgswemv2
```
The docker container contains all dependencies, and builds all targets.

For tips, on how to get started, we recommend looking at our users guide, which can be found in `documentation/users-guide/dgswem-v2-users-guide.pdf`. Basic examples about how to run the code and the design philosophy behind dgswemv2 can be found here.

Lastly, for any further questions, we encourage you to either open an issue on the github repository or join our slack channel. Simply open an issue in the repo with the `Slack` tag and the email address with which you would like to be added.

## Building instructions

### Dependencies

In designing any software package there are clear benefits and disadvantages to using libaries. In dgswem-v2, we have the following dependencies:

| Dependency | Function               |   Dependent Targets       |
| ---------- | ---------------------- | ------------------------- |
| yaml-cpp   | YAML file parser       | All targets               |
| Metis      | Graph partitioner      | `partitioner`             |
| HPX        | Asynchronous Runtime   | `*_SWE_HPX`               |
| OpenMP     | Application threading  | `*_SWE_OMPI`              |
| MPI        | Message passing        | `*_SWE_OMPI`,`*_SWE__HPX` |
| Eigen      | Linear Algebra Library | `EHDG_*`                  |

To begin select a work directory. On the TACC machines we recommend using the `$WORK` directory and copying all of the meshes into `$SCRATCH` for running jobs. `yaml-cpp`, `Metis`, and `HPX` are all libraries. To assist the user in installing the Libraries we have added some bash scripts in `scripts/build` to assist the user in compiling `dgswem-v2`.

To begin the installation process, you will need to set up a configuration file. From the root of the repository, go to `scripts/build`. Open up `config.txt` and define a machine according to your preferences. Note the `build-XXXX.sh` scripts will use `config.txt` by default. However, one can also use custom configuration file as follows:
```sh
    ./build-XXXX.sh -c custom_config.txt
```
Note that for HPX, the library's build-type must be consistent with the applications build-type.

#### Installing `yaml-cpp`

From `$WORK` (as defined in the configuration file)
```sh
    git clone git@github.com:jbeder/yaml-cpp.git
    cd /path/to/build/scripts/
    ./build-yaml-cpp.sh
```
#### Installing `hpx`

The building of hpx requires a lot of dependencies. In particular, for our bash script, we require that jemalloc, boost, hwloc, and an MPI implementation be installed. Note that most clusters typically come installed with boost, hwloc, and an MPI implementation. Thus, if you have already have installed versions of the afore mentioned libraries, skip the relevant build scripts.
```sh
    cd $WORK
    git clone https://github.com/STEllAR-GROUP/hpx
    cd /path/to/build/scripts
    ./build-hwloc.sh
    ./build-boost.sh
    ./build-jemalloc.sh
    ./build-hpx.sh
```

#### Installing `eigen`

To build the hybridized discontinuous Galerkin targets, the application requires Eigen3. Although it is a header only library, we have provided a `build-eigen.sh` script, which will run download Eigen 3.3.4 from the internet, run cmake to run various checks, and lastly, install the headers to the `config.txt` specified location.
```
cd /path/to/build/scripts
./build-eigen.txt
```

### Building the Application

Assuming that `${INSTALL_PATH}` is the path defined in the configuration file. The dgswemv2 uses an out of source build. The application can be built as follows:
```sh
    cd /path/to/dgswemv2/
    mkdir build
    cmake -DCMAKE_BUILD_TYPE=<build type> -DCMAKE_PREFIX_PATH=${INSTALL_PATH} ..
```
Note that there are some additional options, which will create additional targets. These typically require additional dependencies.

| CMake Option   | Description                                                           |
| -------------- | --------------------------------------------------------------------- |
| USE_OMPI       | Enables MPI OpenMP parallelization; builds target `DG_HYPER_SWE_OMPI` |
| USE_HPX        | Enables HPX parallelization; builds target `DG_HYPER_SWE_HPX`         |
| COMPILER_WARNINGS | Display compiler warnings                                          |
| SET_VERBOSE    | Set the makefile compilation output to Verbose                        |
| BUILD_EXAMPLES | Build additional executables to run the examples                      |
| RKDG_SWE       | Build Runge-Kutta discontinuous Galerkin targets                      |
| EHDG_SWE       | Build explicit Hybridized discontinuous Galerkin targets              |
| IHDG_SWE       | Build implicit hybridized discontinuous Galerkin targets              |

Note that by default `RKDG_SWE` is set to `On` and associated targets will be built by cmake.

## License

DGSWEM V2 is licensed under the MIT license. The following files have been copied (and potentially modified) from other repositories. Their licenses are inlined within the files:

 - `cmake/modules/FindMETIS.cmake`
 - `documentation/doxygen-bootstrapped/customdoxygen.css`
 - `documentation/doxygen-bootstrapped/doxy-boot.js`
 - `documentation/doxygen-bootstrapped/footer.html`
 - `documentation/doxygen-bootstrapped/header.html`
 - `source/utilities/linear_algebra/serialization/blaze_vector.hpp`
 - `source/utilities/linear_algebra/serialization/blaze_matrix.hpp`

N.B. The files in `documentation/doxygen-boostrapped` are taken from `feature/support-doxygen-1.1.12+` branch of the [Velron/doxygen-bootstrapped](https://github.com/Velron/doxygen-bootstrapped) repository. The copy of the Apache-2.0 licesnse can be found in `documentation/doxygen-bootstrapped/APACHE_2.0_LICENSE` of this repository or in the original repository.