
add_executable("hi_sdc" FE_sdc.cpp)
target_link_dune_default_libraries("hi_sdc")

add_executable("hi_mlsdc" FE_mlsdc.cpp)
target_link_dune_default_libraries("hi_mlsdc")
