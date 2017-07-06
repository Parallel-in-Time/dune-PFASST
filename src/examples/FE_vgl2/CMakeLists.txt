
add_executable("FE_sdcvgl2" FE_sdcFP.cpp)
target_link_dune_default_libraries("FE_sdcvgl2")

add_executable("FE_mlsdcvgl2" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("FE_mlsdcvgl2")

#add_executable("FE_pfasstvgl2" FE_pfasstFP.cpp)
#target_link_dune_default_libraries("FE_pfasstvgl2")
