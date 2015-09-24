# ===========================================================================
#                  Nexus-Tools
# ===========================================================================

# ----------------------------------------------------------------------------
# Load SeqAn build system
# ----------------------------------------------------------------------------

project (NexusTools)
cmake_minimum_required (VERSION 2.8.2)

if(LINUX_STATIC)
	message("building with linux static libs")
	SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
	SET(BUILD_SHARED_LIBRARIES OFF)
	SET(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} " -static -static-libgcc -static-libstdc++")
    set(Boost_USE_STATIC_LIBS ON)
#    set(LINK_SEARCH_START_STATIC ON)
    set(LINK_SEARCH_END_STATIC ON)
    set(SEQAN_LIBRARIES ${SEQAN_LIBRARIES} "-Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
#	SET_TARGET_PROPERTIES(nexus-pre PROPERTIES LINK_SEARCH_END_STATIC 1)
else()
#	add_definitions (-DBOOST_ALL_DYN_LINK)
endif(LINUX_STATIC)

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

# include(include/seqan/CMakeLists.txt)
seqan_build_system_init ()
seqan_get_repository_info ()
seqan_setup_library ()
include (package)
include (SeqAnCtdSetup)

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

add_definitions (-D_MULTITHREADED_IO)

add_subdirectory(nexus-pre)
add_subdirectory(test-app)
add_subdirectory(ting)
add_subdirectory(MappingAnalyzer)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../seqan/apps/seqan_flexbar ${CMAKE_CURRENT_BINARY_DIR}/flexbar)

if(LINUX_STATIC)
	SET_TARGET_PROPERTIES(nexus-pre PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(test-app PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(ting PROPERTIES LINK_SEARCH_END_STATIC 1)
	SET_TARGET_PROPERTIES(MappingAnalyzer PROPERTIES LINK_SEARCH_END_STATIC 1)
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
