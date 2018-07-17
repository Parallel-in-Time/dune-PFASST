
add_executable("error" error_gridview.cpp)
target_link_dune_default_libraries("error")

add_executable("_sdc" FE_sdcFP.cpp)
target_link_dune_default_libraries("_sdc")

add_executable("_mlsdc" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("_mlsdc")

