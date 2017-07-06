
add_executable("FE_sdcvgl" FE_sdcFP.cpp)
target_link_dune_default_libraries("FE_sdcvgl")

add_executable("FE_mlsdcvgl" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("FE_mlsdcvgl")

add_executable("FE_pfasstvgl" FE_pfasstFP.cpp)
target_link_dune_default_libraries("FE_pfasstvgl")
