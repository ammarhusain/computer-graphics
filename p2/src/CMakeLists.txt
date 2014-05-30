cmake_minimum_required(VERSION 2.8)
project(P2)

include(build/CMakeLists.txt)

add_subdirectory(application)
add_subdirectory(math)
add_subdirectory(p2)
add_subdirectory(scene)

if(APPLE)
	add_subdirectory(SDLmain)
endif()
