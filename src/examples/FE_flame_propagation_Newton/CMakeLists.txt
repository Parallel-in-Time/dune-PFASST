
add_executable("FE_sdcFPN" FE_sdcFP.cpp)
target_link_dune_default_libraries("FE_sdcFPN")

add_executable("FE_mlsdcFPN" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("FE_mlsdcFPN")

#add_executable("FE_pfasstFPN" FE_pfasstFP.cpp)
#target_link_dune_default_libraries("FE_pfasstFPN")
