include_directories(${hybrid_fem_bie_SOURCE_DIR}/src)
add_executable(fish_bone fish_bone.cc)
target_link_libraries(fish_bone hybrid_fem_bie)
