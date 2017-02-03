add_executable("vanilla_sdc" vanilla_sdc.cpp)
target_link_dune_default_libraries("vanilla_sdc")

add_executable("serial_mlsdc" serial_mlsdc.cpp)
target_link_dune_default_libraries("serial_mlsdc")
