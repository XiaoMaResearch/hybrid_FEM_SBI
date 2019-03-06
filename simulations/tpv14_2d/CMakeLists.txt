include_directories(${hybrid_fem_bie_SOURCE_DIR}/src)


add_executable(tpv14_2d_run tpv14_2d_run.cc)
target_link_libraries(tpv14_2d_run hybrid_fem_bie)


add_executable(tpv14_2d_run_with_buffer tpv14_2d_run_with_buffer.cc)
target_link_libraries(tpv14_2d_run_with_buffer hybrid_fem_bie)
