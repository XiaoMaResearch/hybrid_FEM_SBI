include_directories(${hybrid_fem_bie_SOURCE_DIR}/src)



add_executable(plasticity_test parallel_fault_run_reg_plastic.cc)
target_link_libraries(plasticity_test hybrid_fem_bie)
