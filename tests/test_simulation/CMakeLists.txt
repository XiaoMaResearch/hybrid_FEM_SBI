include_directories(${hybrid_fem_bie_SOURCE_DIR}/src)

# Multi faults system (Trying to get it general)



add_executable(input_Slip_Weakening_Generated_multi_faults input_Slip_Weakening_Generated_multi_faults.cc)
target_link_libraries(input_Slip_Weakening_Generated_multi_faults hybrid_fem_bie)
