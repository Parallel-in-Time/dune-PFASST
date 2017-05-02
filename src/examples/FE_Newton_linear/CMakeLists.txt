
add_executable("FE_sdcNFP" FE_sdcFP.cpp)
target_link_dune_default_libraries("FE_sdcNFP")

add_executable("FE_mlsdcNFP" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("FE_mlsdcNFP")

#add_executable("FE_pfasstNFP" FE_pfasstFP.cpp)
#target_link_dune_default_libraries("FE_pfasstNFP")
