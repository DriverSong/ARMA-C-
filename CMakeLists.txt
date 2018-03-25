cmake_minimum_required(VERSION 3.9)
project(MyARMA)

set(CMAKE_CXX_STANDARD 11)

add_executable(MyARMA main.cpp ARMACore.cpp ARMACore.h AR.h MA.h ARMA.h Predict.cpp Predict.h Control.cpp)