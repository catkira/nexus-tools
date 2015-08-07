cmake_minimum_required(VERSION 3.1)
project(nexus-pre)
message (STATUS "Configuring nexus-pre")

set (SEQAN_FIND_DEPENDENCIES OpenMP ZLIB BZip2)
find_package(SeqAn REQUIRED)
find_package(CXX11 REQUIRED)
find_package(CXX14)

if (NOT CXX11_FOUND)
    message (FATAL_ERROR "  C++11 not found: Not building lambda.")
    return ()
endif (NOT CXX11_FOUND)

# Warn if OpenMP was not found.
if (NOT SEQAN_HAS_OPENMP)
    message (WARNING "WARNING WARNING WARNING\nWARNING: OpenMP not found. nexus-pre will be built without multi-threading! "
    "This is probably not what you want! Use GCC or a very, very recent version of Clang.\nWARNING WARNING WARNING")
endif (NOT SEQAN_HAS_OPENMP)

# Warn if Zlib was not found.
if (NOT SEQAN_HAS_ZLIB)
  message (WARNING "WARNING: Zlib not found. Building nexus-pre without gzipped input and output.")
endif (NOT SEQAN_HAS_ZLIB)


# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Sources files
    FILE(GLOB Sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp) 
	
# Header files
    FILE(GLOB Headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h) 
	

#include(${Projects_SOURCE_DIR}/SourceGroups.cmake)	
		
# Output
    SET(ALLSRC ${Sources} ${Headers})
    ADD_EXECUTABLE(nexus-pre ${ALLSRC})
	
    target_link_libraries (nexus-pre ${SEQAN_LIBRARIES})
				  
