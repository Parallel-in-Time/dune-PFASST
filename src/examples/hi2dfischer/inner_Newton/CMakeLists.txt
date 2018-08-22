
add_executable("sdc" sdc.cpp)
target_link_dune_default_libraries("sdc")

add_executable("mlsdc" mlsdc.cpp)
target_link_dune_default_libraries("mlsdc")

add_executable("pfasst" pfasst.cpp)
target_link_dune_default_libraries("pfasst")
