cmake_minimum_required(VERSION 3.17)
project(vmlab6)

set(CMAKE_CXX_STANDARD 20)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

set(BOOST_INCLUDE_DIRS "/opt/local/include")
set(BOOST_LIBRARY_DIRS "/opt/local/lib")

FIND_PACKAGE(Boost REQUIRED COMPONENTS filesystem iostreams)

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
add_executable(vmlab6 main.cpp ODE.cpp ODE.h)
target_link_libraries(vmlab6 ${Boost_LIBRARIES})