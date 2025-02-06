cmake_minimum_required(VERSION 3.24)
project(example-cmd VERSION 1.0.0)
set(CMAKE_CXX_STANDARD 11)
add_executable(stretch ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

include(FetchContent)
set(FETCHCONTENT_QUIET False)
FetchContent_Declare(
	inputs
	URL https://signalsmith-audio.co.uk/code/stretch/inputs.zip
	SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/inputs"
)
FetchContent_MakeAvailable(inputs)
