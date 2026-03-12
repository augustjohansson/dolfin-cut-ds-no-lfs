cmake_minimum_required(VERSION 3.5)
project(dolfin-tests)

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

# Prepare "Catch" library for other executables
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/catch)

add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

# New geometry tests (use relative "../geometry/..." includes and local
# "catch/catch.hpp").  The include directory ${DOLFIN_SOURCE_DIR}/dolfin/geometry
# is added below so that "../geometry/Foo.h" resolves to the DOLFIN geometry
# source directory.
set(GEOMETRY_TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/geometry/testCollisionPredicates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/geometry/testConvexTriangulation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/geometry/testIntersectionConstruction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/geometry/testSimplexQuadrature.cpp
  # testCGALComparison.cpp is guarded by #ifdef DOLFIN_WITH_CGAL internally;
  # always include it in the build – it compiles to an empty TU when CGAL is off.
  ${CMAKE_CURRENT_SOURCE_DIR}/geometry/testCGALComparison.cpp
  )

# Make test executable
set(TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/common/SubSystemsManager.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/function/Expression.cpp
  ${GEOMETRY_TEST_SOURCES}
  ${CMAKE_CURRENT_SOURCE_DIR}/io/XMLMeshData.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/io/XMLMeshValueCollection.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/la/LinearOperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/la/Vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/Mesh.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/MeshColoring.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/MeshFunction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/MeshValueCollection.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/parameter/Parameters.cpp
  )

add_executable(unittests ${TEST_SOURCES})
target_link_libraries(unittests PRIVATE Catch dolfin)

# Primary include directories: DOLFIN source tree and build tree.
# Also add ${DOLFIN_SOURCE_DIR}/dolfin/geometry so that the relative
# include path "../geometry/Foo.h" used in the new geometry test files
# resolves correctly to ${DOLFIN_SOURCE_DIR}/dolfin/geometry/Foo.h.
target_include_directories(unittests PRIVATE
  ${DOLFIN_SOURCE_DIR}
  ${DOLFIN_SOURCE_DIR}/dolfin
  ${DOLFIN_SOURCE_DIR}/dolfin/geometry
  ${DOLFIN_BINARY_DIR}
  )

# When CGAL is enabled the geometry test files indirectly include CGAL
# headers (via CGALPredicates.h).  Add CGAL include directories so the
# compiler can find them.  Also link against GMP and MPFR which are
# required at link time by CGAL's exact-constructions kernel (EPECK) used
# in testCGALComparison.cpp; libdolfin links them PRIVATE so they are not
# transitively available to dependents.
if (DOLFIN_USE_CGAL AND CGAL_FOUND)
  target_include_directories(unittests PRIVATE ${CGAL_INCLUDE_DIRS})
  if (GMP_LIBRARIES)
    target_link_libraries(unittests PRIVATE ${GMP_LIBRARIES})
  endif()
  if (MPFR_LIBRARIES)
    target_link_libraries(unittests PRIVATE ${MPFR_LIBRARIES})
  endif()
endif()

# Test target
add_test(unittests unittests)
