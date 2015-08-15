cmake_minimum_required(VERSION 3.1)
project(nexus-pre)
message (STATUS "Configuring nexus-pre")



# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------




# Sources files
    FILE(GLOB Sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp) 
	
# Header files
    FILE(GLOB Headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h) 
	

#include(${Projects_SOURCE_DIR}/SourceGroups.cmake)	
		
# Output
    SET(ALLSRC ${Sources} ${Headers})
    ADD_EXECUTABLE(nexus-pre ${ALLSRC})
	
    target_link_libraries (nexus-pre ${SEQAN_LIBRARIES})
				  
