Stomi
===================

Stomi is a toolbox that can be used to study stochastic migration procesess in the Solar System. Stomi models the process of stochastic migration by simplifying the gravitational N-body problem as a superposition of multiple 3-body problems followed by a random walk.

This toolbox can be used to study the stochastic migration of a massive body, orbiting a central body, due to the gravitational interactions with a ring of low-mass perturbers. 

Possible applications include:
- Studying migration of Neptune due to interactions with the Kuiper belt
- Studying migration of Ceres and Vesta due to interactions with the Main Asteroid belt
- Studying migration of late-stage protoplanets due to scattering interactions with minor bodies
- Studying migration of exoplanets embedded in a ring of low-mass bodies 

Requirements
------

`stomi` requires the following libraries/tools. The listed versions are not necessarily essential, however `Stomi` is known to build succesfully when these are used.

| Name                                                                     | Version       |
| -------------                                                            |:-------------:|
| [GCC](http://gcc.gnu.org "GCC homepage")                                 | 4.8           |
| [CMake](http://www.cmake.org/ "CMake homepage")                        | 2.8.12        |
| [Eigen](http://eigen.tuxfamily.org "Eigen's homepage")                 | 3.2.1         |
| [Boost](http://www.boost.org "Boost's homepage")                       | 1.55.0        |
| [Tudat & Tudat Core](http://tudat.tudelft.nl "Tudat project homepage") | [2a2720](https://github.com/kartikkumar/tudat-svn-mirror/tree/54dc69cd91e84c2a9cddc4caf9f0e86aba2a2720) & [bec885](https://github.com/kartikkumar/tudatCore-svn-mirror/tree/184a180d7213aeb021d672b7b92b0733a4bec885) |
| [Assist](https://github.com/kartikkumar/assist "Assist project")                       | [f9c915](https://github.com/kartikkumar/assist/tree/c3f8281dc21d0d7364aecd63c8ea68e929f9c915)        |

Alternative C++ compilers may be used, however they haven't been tested to date.

Installation
------


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
