
add_executable("par" parallelCG.cpp)
target_link_dune_default_libraries("par")

add_executable("sdc" FE_sdcFP.cpp)
target_link_dune_default_libraries("sdc")

add_executable("mlsdc" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("mlsdc")

add_executable("pfasst" FE_pfasstFP.cpp)
target_link_dune_default_libraries("pfasst")
