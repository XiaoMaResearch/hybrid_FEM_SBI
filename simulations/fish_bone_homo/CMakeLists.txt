include_directories(${hybrid_fem_bie_SOURCE_DIR}/src)

add_executable(fish_bone_homo fish_bone_homo.cc)
target_link_libraries(fish_bone_homo hybrid_fem_bie)
