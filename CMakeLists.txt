# ===========================================================================
#                  Lambda
# ===========================================================================

# ----------------------------------------------------------------------------
# Load SeqAn build system
# ----------------------------------------------------------------------------

project (seqan)
cmake_minimum_required (VERSION 2.8.2)
set (CMAKE_INCLUDE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/../seqan/include
    ${CMAKE_INCLUDE_PATH})
set (CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/../seqan/util/cmake/FindTBB
    ${CMAKE_CURRENT_SOURCE_DIR}/../seqan/util/cmake
    ${CMAKE_MODULE_PATH})
include (SeqAnContribs)
set (SEQAN_USE_SEQAN_BUILD_SYSTEM TRUE CACHE INTERNAL "Use SeqAn build system." FORCE)
include (SeqAnBuildSystem)
set (SEQAN_ROOT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}../seqan" CACHE INTERNAL "Root source directory." FORCE)
option(LINUX_STATIC "build with static libs" OFF)

if(LINUX_STATIC)
	message("building with linux static libs")
	SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
	SET(BUILD_SHARED_LIBRARIES OFF)
	SET(CMAKE_EXE_LINKER_FLAGS "-static")
#	SET_TARGET_PROPERTIES(nexus-pre PROPERTIES LINK_SEARCH_END_STATIC 1)
endif(LINUX_STATIC)

# include(include/seqan/CMakeLists.txt)
seqan_build_system_init ()
seqan_get_repository_info ()
seqan_setup_library ()
include (package)
include (SeqAnCtdSetup)
find_package (SeqAn REQUIRED)
find_package (CXX11 REQUIRED)
find_package (CXX14 REQUIRED)

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS} ${CXX14_CXX_FLAGS}")

add_subdirectory(nexus-pre)
add_subdirectory(test-app)
add_subdirectory(ting)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../seqan/apps/seqan_flexbar ${CMAKE_CURRENT_BINARY_DIR}/flexbar)

if(LINUX_STATIC)
	SET_TARGET_PROPERTIES(nexus-pre PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(test-app PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(ting PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(seqan_flexlib PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(seqan_flexbar PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(sflexAR PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(sflexDMulti PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(sflexFilter PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(sflexQC PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(test_sflex_demultiplex PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(test_sflex_general_processing PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(test_sflex_trimming PROPERTIES LINK_SEARCH_END_STATIC 1)
endif(LINUX_STATIC)

