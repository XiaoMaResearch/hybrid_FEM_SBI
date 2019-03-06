include_directories(${hybrid_fem_bie_SOURCE_DIR}/src)


add_executable(tpv5_2d tpv5_2d.cc)
target_link_libraries(tpv5_2d hybrid_fem_bie)

# add_executable(tpv5_2d_run tpv5_2d_run.cc)
# target_link_libraries(tpv5_2d_run hybrid_fem_bie)
#
# add_executable(tpv5_2d_generated_run tpv5_2d_generated_run.cc)
# target_link_libraries(tpv5_2d_generated_run hybrid_fem_bie)
