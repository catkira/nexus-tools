cmake_minimum_required(VERSION 3.1)
project(test-app)
message (STATUS "Configuring test-app")

add_definitions (-DSEQAN_APP_VERSION="1.0.0")

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES OpenMP ZLIB)
find_package (SeqAn REQUIRED)
find_package (CXX11 REQUIRED)
find_package (CXX14 REQUIRED)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS} ${CXX14_CXX_FLAGS}")

# Warn if OpenMP was not found.
if (NOT ${SEQAN_HAS_OPENMP})
    message (WARNING "WARNING WARNING WARNING\nWARNING: OpenMP not found. test-app will be built without multi-threading! "
    "This is probably not what you want! Use GCC or a very, very recent version of Clang.\nWARNING WARNING WARNING")
endif (NOT ${SEQAN_HAS_OPENMP})

# Warn if Zlib was not found.
if (NOT ${SEQAN_HAS_ZLIB})
  message (WARNING "WARNING: Zlib not found. Building test-app without gzipped input and output.")
endif (NOT ${SEQAN_HAS_ZLIB})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------
SET(Boost_USE_STATIC_LIBS ON)
FIND_PACKAGE( Boost 1.40 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES(${SEQAN_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../shared_headers)
LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS} )

# Sources files
    FILE(GLOB Sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp) 
	
# Header files
    FILE(GLOB Headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h) 
	set(Headers ${Headers} ${CMAKE_CURRENT_SOURCE_DIR}/../shared_headers/BamRecordKey.h)

#include(${Projects_SOURCE_DIR}/SourceGroups.cmake)	
		
# Output
    SET(ALLSRC ${Sources} ${Headers})
    ADD_EXECUTABLE(test-app ${ALLSRC})
	
    target_link_libraries (test-app ${SEQAN_LIBRARIES} ${Boost_LIBRARIES} )
				  
