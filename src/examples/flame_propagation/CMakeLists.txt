
add_executable("sdcFP" FE_sdcFP.cpp)
target_link_dune_default_libraries("sdcFP")

add_executable("mlsdcFP" FE_mlsdcFP.cpp)
target_link_dune_default_libraries("mlsdcFP")

#add_executable("pfasstFP" FE_pfasstFP.cpp)
#target_link_dune_default_libraries("pfasstFP")
