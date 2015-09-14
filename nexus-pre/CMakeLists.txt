cmake_minimum_required(VERSION 3.1)
project(nexus-pre)
message (STATUS "Configuring nexus-pre")

add_definitions (-DSEQAN_APP_VERSION="1.0.0")

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES OpenMP ZLIB)
find_package (SeqAn REQUIRED)

# Warn if OpenMP was not found.
if (NOT ${SEQAN_HAS_OPENMP})
    message (WARNING "WARNING WARNING WARNING\nWARNING: OpenMP not found. nexus-pre will be built without multi-threading! "
    "This is probably not what you want! Use GCC or a very, very recent version of Clang.\nWARNING WARNING WARNING")
endif (NOT ${SEQAN_HAS_OPENMP})

# Warn if Zlib was not found.
if (NOT ${SEQAN_HAS_ZLIB})
  message (WARNING "WARNING: Zlib not found. Building nexus-pre without gzipped input and output.")
endif (NOT ${SEQAN_HAS_ZLIB})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------
FIND_PACKAGE( Boost 1.40 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES(${SEQAN_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../shared_headers)



# Sources files
    FILE(GLOB Sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp) 
	
# Header files
    FILE(GLOB Headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h) 
	set(Headers ${Headers} ${CMAKE_CURRENT_SOURCE_DIR}/../shared_headers/BamRecordKey.h)
	set(Headers ${Headers} ${CMAKE_CURRENT_SOURCE_DIR}/../shared_headers/peak.h)

#include(${Projects_SOURCE_DIR}/SourceGroups.cmake)	
		
# Output
    SET(ALLSRC ${Sources} ${Headers})
    ADD_EXECUTABLE(nexus-pre ${ALLSRC})
	
    target_link_libraries (nexus-pre ${SEQAN_LIBRARIES} ${Boost_LIBRARIES} )
				  
