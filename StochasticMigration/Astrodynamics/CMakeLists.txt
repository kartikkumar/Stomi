 #
 #    Copyright (c) 2010-2013, Delft University of Technology
 #    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 #    All rights reserved.
 #    See http://bit.ly/12SHPLR for license details.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      120820    K. Kumar          Adapted file to work with MabSimulations/Astrodynamics.
 #      130217    K. Kumar          Updated static library prefix to stochastic_migration_.
 #
 #    References
 #
 #    Notes
 #

# Set source files.
set(ASTRODYNAMICS_SOURCES
  "${SRCROOT}${ASTRODYNAMICSDIR}/propagationDataPoint.cpp"
  "${SRCROOT}${ASTRODYNAMICSDIR}/randomWalkFunctions.cpp"
)

# Set header files.
set(ASTRODYNAMICS_HEADERS
  "${SRCROOT}${ASTRODYNAMICSDIR}/hillSphere.h"
  "${SRCROOT}${ASTRODYNAMICSDIR}/propagationDataPoint.h"
  "${SRCROOT}${ASTRODYNAMICSDIR}/randomWalkFunctions.h"
)

# Add static library.
add_library(stochastic_migration_astrodynamics STATIC ${ASTRODYNAMICS_SOURCES} ${ASTRODYNAMICS_HEADERS})
setup_library_target(stochastic_migration_astrodynamics "${SRCROOT}${ASTRODYNAMICSDIR}")
