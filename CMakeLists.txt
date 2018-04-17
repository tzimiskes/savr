cmake_minimum_required(VERSION 2.8)
project(savr)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(R REQUIRED)
message("R command: ${R_COMMAND}")
message("Rscript command: ${RSCRIPT_COMMAND}")
message("R home: ${R_HOME}")
message("R libraries: ${R_LIBRARIES}")
message("R include dirs: ${R_INCLUDE_DIR} ${Rcpp_INCLUDE_DIR}")

find_package(savvy REQUIRED)

set(CMAKE_CXX_STANDARD 11)

add_library(savr SHARED src/savr.cpp)
include_directories(savr ${R_INCLUDE_DIR} ${Rcpp_INCLUDE_DIR})
target_link_libraries(savr ${R_LIBRARIES} savvy)

install(TARGETS savr DESTINATION lib)
