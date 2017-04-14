add_executable("sdcGS" FE_sdcGS.cpp)
target_link_dune_default_libraries("sdcGS")

add_executable("mlsdcGS" FE_mlsdcGS.cpp)
target_link_dune_default_libraries("mlsdcGS")

add_executable("pfasstGS" FE_pfasstGS.cpp)
target_link_dune_default_libraries("pfasstGS")
