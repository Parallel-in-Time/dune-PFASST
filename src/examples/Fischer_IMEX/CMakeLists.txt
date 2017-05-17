
add_executable("sdc_imex" FE_sdcFP.cpp)
target_link_dune_default_libraries("sdc_imex")

add_executable("mlsdc_imex" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("mlsdc_imex")

#add_executable("FE_pfasstNFP" FE_pfasstFP.cpp)
#target_link_dune_default_libraries("FE_pfasstNFP")
