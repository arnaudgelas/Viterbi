project( Viberti )

cmake_minimum_required( VERSION 2.8 )

find_package( Boost REQUIRED )

include_directories( ${Boost_INCLUDE_DIRS} )
include(CTest)

add_executable( ViterbiTest test.cpp )
add_test( NAME ViterbiTest COMMAND ${CMAKE_BINARY_DIR}/ViterbiTest )

add_executable( ViterbiTest2 test2.cpp )
add_test( NAME ViterbiTest2 COMMAND ${CMAKE_BINARY_DIR}/ViterbiTest2 )
