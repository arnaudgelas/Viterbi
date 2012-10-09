project( Viberti )

cmake_minimum_required( VERSION 2.8 )

find_package( Boost REQUIRED )

include_directories( ${Boost_INCLUDE_DIRS} )

add_executable( viterbi test.cpp )
