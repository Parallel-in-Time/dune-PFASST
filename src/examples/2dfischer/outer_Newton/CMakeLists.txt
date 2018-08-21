add_executable("vor" vor_pfasstFP.cpp)
target_link_dune_default_libraries("vor")

add_executable("newton_sdc" FE_sdcFP.cpp)
target_link_dune_default_libraries("newton_sdc")

add_executable("newton_mlsdc" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("newton_mlsdc")

add_executable("newton_pfasst" FE_pfasstFP.cpp)
target_link_dune_default_libraries("newton_pfasst")


