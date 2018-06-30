
add_executable("run" run.cpp)
target_link_dune_default_libraries("run")

add_executable("_sdc" FE_sdcFP.cpp)
target_link_dune_default_libraries("_sdc")

add_executable("_mlsdc" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("_mlsdc")

add_executable("_pfasst" FE_pfasstFP.cpp)
target_link_dune_default_libraries("_pfasst")
