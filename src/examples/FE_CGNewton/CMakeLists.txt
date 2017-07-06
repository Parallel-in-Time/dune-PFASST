
add_executable("FE_sdcCG" FE_sdcFP.cpp)
target_link_dune_default_libraries("FE_sdcCG")

add_executable("FE_mlsdcCG" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("FE_mlsdcCG")


