cmake_minimum_required(VERSION 2.8)

if(MSVC)
	add_definitions(-D__SIZEOF_FLOAT__=4 -D__SIZEOF_FLOAT__=4 -D__SIZEOF_DOUBLE__=8 -Drestrict=__restrict -D_USE_MATH_DEFINES)
endif(MSVC)

SET(SOURCE current.c emf.c particles.c random.c timer.c main.c simulation.c zdf.c)
add_executable(zpic ${SOURCE})
