include_directories(${hybrid_fem_bie_SOURCE_DIR}/src)
link_directories(~/Cpp_Libarary/seacas/lib)
#add_executable(biegel_test biegel_test.cc)
#target_link_libraries(biegel_test hybrid_fem_bie)

add_executable(biegel_test_plastic_exodus biegel_test_plastic_exodus.cc)
target_link_libraries(biegel_test_plastic_exodus hybrid_fem_bie libexodus.a libnetcdf.dylib)
