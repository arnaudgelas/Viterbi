project( Viberti )

cmake_minimum_required( VERSION 2.8 )

find_package( Boost REQUIRED )

include_directories( ${Boost_INCLUDE_DIRS} )

add_executable( ViterbiTest test.cpp )
add_executable( ViterbiTest2 test2.cpp )
