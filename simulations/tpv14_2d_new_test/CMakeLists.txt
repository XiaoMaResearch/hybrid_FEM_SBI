include_directories(${hybrid_fem_bie_SOURCE_DIR}/src)


add_executable(tpv14_2d_new_test tpv14_2d_new_test.cc)
target_link_libraries(tpv14_2d_new_test hybrid_fem_bie)
