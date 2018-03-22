cmake_minimum_required(VERSION 3.9)
project(MyARMA)

set(CMAKE_CXX_STANDARD 11)

add_executable(MyARMA main.cpp ARMACore.cpp ARMACore.h AR.cpp AR.h MA.cpp MA.h ARMA.cpp ARMA.h)