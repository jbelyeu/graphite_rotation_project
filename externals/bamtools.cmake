#  For more information, please see: http://software.sci.utah.edu
#
#  The MIT License
#
#  Copyright (c) 2015 Scientific Computing and Imaging Institute,
#  University of Utah.
#
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.

# Setting up external library bamtools, we don't build it because we only need the include directories

SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(BAMTOOLS_PROJECT bamtools_project CACHE INTERNAL "bamtools project name")
SET(BAMTOOLS_DIR ${CMAKE_BINARY_DIR}/externals/bamtools CACHE INTERNAL "bamtools project directory")
SET(BAMTOOLS_LIB)
ExternalProject_Add(${BAMTOOLS_PROJECT}
	GIT_REPOSITORY https://github.com/pezmaster31/bamtools.git
	GIT_TAG 2d7685d2aeedd11c46ad3bd67886d9ed65c30f3e #lock in the commit id so we don't this doesn't break in the future
	INSTALL_COMMAND ""
    UPDATE_COMMAND ""
	PREFIX ${BAMTOOLS_DIR}
    CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
	CMAKE_ARGS
		-DCMAKE_BUILD_TYPE=Debug
)

ExternalProject_Get_Property(${BAMTOOLS_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${BAMTOOLS_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${BAMTOOLS_PROJECT} BINARY_DIR)

SET(BAMTOOLS_LIB ${SOURCE_DIR}/lib/libbamtools.a CACHE INTERNAL "Bamtools Lib")
SET(BAMTOOLS_UTIL_LIB ${SOURCE_DIR}/lib/libbamtools-utils.a CACHE INTERNAL "Bamtools-util Lib")
SET(BAMTOOLS_INCLUDE ${SOURCE_DIR}/include CACHE INTERNAL "Bamtools Include")
