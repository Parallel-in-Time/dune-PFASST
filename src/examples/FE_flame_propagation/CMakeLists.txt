
add_executable("FE_sdcFP" FE_sdcFP.cpp)
target_link_dune_default_libraries("FE_sdcFP")

add_executable("FE_mlsdcFP" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("FE_mlsdcFP")

#add_executable("FE_pfasstFP" FE_pfasstFP.cpp)
#target_link_dune_default_libraries("FE_pfasstFP")
