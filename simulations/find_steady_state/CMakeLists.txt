include_directories(${hybrid_fem_bie_SOURCE_DIR}/src)
link_directories(~/Cpp_Libarary/seacas/lib)

add_executable(homo homo.cc)
target_link_libraries(homo hybrid_fem_bie libexodus.a libnetcdf.dylib)
