 #
 #    Copyright (c) 2010-2013, Delft University of Technology
 #    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 #    All rights reserved.
 #    See http://bit.ly/12SHPLR for license details.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      130212    K. Kumar          Updated to work with Basics sub-directory.
 #      130217    K. Kumar          Updated static library prefix to stochastic_migration_.
 #
 #    References
 #
 #    Notes
 #

# Set source files.
set(BASICS_SOURCES
  "${SRCROOT}${BASICSDIR}/basics.cpp"
)

# Set header files.
set(BASICS_HEADERS
  "${SRCROOT}${BASICSDIR}/basics.h"
)

# Add static library.
add_library(stochastic_migration_basics STATIC ${BASICS_SOURCES} ${BASICS_HEADERS})
setup_library_target(stochastic_migration_basics "${SRCROOT}${BASICSDIR}")
