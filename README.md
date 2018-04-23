# BES

Basic Euler Solver, a lightweight, multi-threaded, three-dimensional unstructured finite volume code to solve the Euler equations.

BES uses OpenMP, Eigen and AMGCL.

## Getting Started

Before deciding to compile the code, you can check out the binaries for [Windows 64bit](https://drive.google.com/file/d/1NRtUkby0Tfmqq1eHPAVDsC-AbSJu2N-I/view?usp=sharing) or [Linux 64 bit](). 

If you want to compile, you just need to have the source file for Boost somewhere on your machine. Some header-only libraries from Boost is used in AMGCL, hence the dependency. The required Eigen and AMGCL libraries are included in BES include directory, so no action is needed for that.

### Prerequisites

You can download Boost [here](https://www.boost.org/users/history/). Boost v1.66.0 was used in the development of BES.


### Installing

The compilation of BES is from a simple makefile, Automake or Autoconf was not used. As such, you will need to manually input some options into the makefile.

First, enter your OS in line 2 of the makefile. For Windows:

```
OS = WINDOWS
```

For Linux:

```
OS = LINUX
```

Second, enter the path to Boost in your machine in line 36 of the makefile:

```
BOOST = path/to/boost
```

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc