project(abm LANGUAGES CXX)

# Require C++14-compliant compiler; only available for CMake v. 3.1 and up
set(CMAKE_CXX_STANDARD 14)

cmake_minimum_required(VERSION 3.1)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# General compile settings
IF (NOT CMAKE_BUILD_TYPE)
 SET(CMAKE_BUILD_TYPE "Debug")
 #SET(CMAKE_BUILD_TYPE "Release")
ENDIF (NOT CMAKE_BUILD_TYPE)

# GNU specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()

# Intel specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

# Clang specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
endif()

# CMake seems to have no way to enable/disable testing per subproject,
# so we provide an option similar to BUILD_TESTING, but just for ABM.
option(ABM_BUILD_TESTING "enable testing for abm" ON)

# CMake Modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# MPI
find_package(MPI)
if (MPI_FOUND)
  add_definitions("-DUSE_MPI")
  include_directories(${MPI_CXX_INCLUDE_DIRS})
  link_libraries(${MPI_CXX_LIBRARIES})
endif()

# OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Include directories
include_directories(BEFORE
  ${abm_SOURCE_DIR}/include/
  ${abm_SOURCE_DIR}/external/
)

# abm executable
SET(abm_src
  ${abm_SOURCE_DIR}/src/graph.cc
  ${abm_SOURCE_DIR}/src/router.cc
  ${abm_SOURCE_DIR}/src/agent.cc
  ${abm_SOURCE_DIR}/src/router_hybrid.cc
)

add_executable(abm ${abm_src} ${abm_SOURCE_DIR}/src/main.cc)

# Unit test
if(ABM_BUILD_TESTING)
  SET(test_src
    ${abm_SOURCE_DIR}/tests/test_main.cc
    ${abm_SOURCE_DIR}/tests/graph_test.cc
    ${abm_SOURCE_DIR}/tests/router_test.cc
    ${abm_SOURCE_DIR}/tests/agent_test.cc
    )
  add_executable(abmtest ${abm_src} ${test_src})
  add_test(NAME abmtest COMMAND $<TARGET_FILE:abmtest>)
  enable_testing()
endif()

# Coverage
find_package(codecov)
if(ENABLE_COVERAGE)
  #add_executable(abmtest_coverage ${abm_src} ${test_src})
  add_coverage(abmtest)
endif()
