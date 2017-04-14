
add_executable("FE_sdcGF" FE_sdcFP.cpp)
target_link_dune_default_libraries("FE_sdcGF")

add_executable("FE_mlsdcGF" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("FE_mlsdcGF")

#add_executable("FE_pfasstFP" FE_pfasstFP.cpp)
#target_link_dune_default_libraries("FE_pfasstFP")
