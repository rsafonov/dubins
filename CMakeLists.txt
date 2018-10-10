cmake_minimum_required(VERSION 2.8)
PROJECT(dubins)
SET(CMAKE_BUILD_TYPE "Debug")

include_directories(include ../dubins/include)

add_library(dubinscar SHARED
  src/dubinscar.cpp
)

set(PKG_NAME "dubins")
set(PKG_DESC "Dubins Library")
set(PKG_DUBINS_LIBS "-ldubinscar")
set(pkg_conf_file "${CMAKE_CURRENT_SOURCE_DIR}/dubins.pc")
configure_file("${pkg_conf_file}.in" "${pkg_conf_file}" @ONLY)

install(FILES ${pkg_conf_file} DESTINATION lib/pkgconfig/ COMPONENT pkgconfig)

install(DIRECTORY include/dubins DESTINATION include
	FILES_MATCHING PATTERN "*.h")

install(TARGETS dubinscar
	LIBRARY DESTINATION lib)

add_executable(dubinstest src/dubinstest.cpp)
target_link_libraries(dubinstest dubinscar tinyxml)


