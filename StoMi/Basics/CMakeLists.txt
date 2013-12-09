 #    Copyright (c) 2010-2013, Delft University of Technology
 #    Copyright (c) 2010-2013, K. Kumar (me@kartikkumar.com)
 #    All rights reserved.
 #    See http://bit.ly/12SHPLR for license details.

# Set source files.
set(BASICS_SOURCES
  "${SRCROOT}${BASICSDIR}/basics.cpp"
)

# Set header files.
set(BASICS_HEADERS
  "${SRCROOT}${BASICSDIR}/basics.h"
)

# Add static library.
add_library(stomi_basics STATIC ${BASICS_SOURCES} ${BASICS_HEADERS})
setup_library_target(stomi_basics "${SRCROOT}${BASICSDIR}")
