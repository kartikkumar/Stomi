Stomi
===================

Stomi is a toolbox that can be used to study stochastic migration procesess in the Solar System. The toolbox breaks down the full N-body problem into a two-stage sequence:

1. Test particle simulator, implementing the Restricted 3-Body Problem
2. Random walk simulator, aggregating the back-reaction of a cloud of perturbers (using the test particle simulation data) on the large migrating body

For an example use of this toolbox, take a look at (Kumar et al., 2014).

Stomi can be used to study the stochastic migration of a massive body, orbiting a central body, due to the gravitational interactions with a ring of low-mass perturbers. 

Possible applications include:
- Studying migration of Neptune due to interactions with the Kuiper belt
- Studying migration of Ceres and Vesta due to interactions with the Main Asteroid belt
- Studying migration of late-stage protoplanets due to scattering interactions with minor bodies
- Studying migration of exoplanets embedded in a ring of low-mass bodies 

Requirements
------

`Stomi` requires the following libraries/tools. The listed versions are not necessarily essential, however `Stomi` is known to build succesfully when these are used.

| Name                                                                     | Version       |
| -------------                                                            |:-------------:|
| [GCC](http://gcc.gnu.org "GCC homepage")                                 | 4.8           |
| [CMake](http://www.cmake.org/ "CMake homepage")                          | 2.8.12        |
| [Eigen](http://eigen.tuxfamily.org "Eigen's homepage")                   | 3.2.1         |
| [Boost](http://www.boost.org "Boost's homepage")                         | 1.55.0        |
| [SQLite](https://sqlite.org/ "SQLite homepage")                          | 3.7.11        |
| [SQLiteCpp](http://srombauts.github.com/SQLiteCpp "SQLite project")      | [d42964](https://github.com/kartikkumar/SQLiteCpp/tree/60b6593998e396f14010e8c82618a7196fd42964)    |
| [Tudat & Tudat Core](http://tudat.tudelft.nl "Tudat project homepage")   | [2a2720](https://github.com/kartikkumar/tudat-svn-mirror/tree/54dc69cd91e84c2a9cddc4caf9f0e86aba2a2720) 
  & [bec885](https://github.com/kartikkumar/tudatCore-svn-mirror/tree/184a180d7213aeb021d672b7b92b0733a4bec885)                                    |
| [Assist](https://github.com/kartikkumar/assist "Assist project")         | [f9c915](https://github.com/kartikkumar/assist/tree/c3f8281dc21d0d7364aecd63c8ea68e929f9c915)       |

Alternative C++ compilers may be used, however they haven't been tested to date.

Installation
------

The easiest way to install Stomi is to simply clone the repository and run `CMake`. This assumes that the above mentioned libraries/tools have been installed and can either be found through the system path or through relative paths (see `CMakeModules` directory included in the Stomi project).

In addition, you must have a C++ compiler and the [CMake](http://www.cmake.org/ "CMake homepage") tool needs to be installed on the system path. 

To clone the `Stomi` toolbox, simply execute the following from your terminal:

```
> git clone https://github.com/kartikkumar/Stomi.git
```

Move into the `Stomi` directory and create a `build` directory:

```
> mkdir build
> cd build
```

From within the `build` directory, run the following command:

```
build> cmake ../Stomi
```

Once `CMake` has been run succesfully, run `make` from within the `build` directory to build the library. 

```
build> make
```

You can optionally also run the test suite at this stage to make sure that everything has been built correctly.

```
build> make test
```

This should not result in any failures.

The applications reside in `bin/applications` directory, relative to the project root.

Documentation
-------------

You can pass the `-DBUILD_DOCUMENTATION=on` option to `CMake` to build the [Doxygen](http://www.doxygen.org "Doxygen homepage") documentation locally.

Contributing
------------

Once you've made your great commits:

1. [Fork](https://github.com/kartikkumar/Stomi/fork) Stomi
2. Create a topic branch - `git checkout -b my_branch`
3. Push to your branch - `git push origin my_branch`
4. Create a [Pull Request](http://help.github.com/pull-requests/) from your
   branch
5. That's it!

License
------

See [COPYING](https://github.com/kartikkumar/Stomi/blob/master/COPYING).

Disclaimer
------

I am not liable for silly use of the Stomi project.

Contact
------

You can reach me at [me@kartikkumar.com](me@kartikkumar.com).
