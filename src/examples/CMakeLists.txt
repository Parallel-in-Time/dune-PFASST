#add_subdirectory(FE_heat)
#add_subdirectory(fischer)
add_subdirectory(Newton)

#add_executable("parallel_solver" parallel_solver.cpp)
#target_link_dune_default_libraries("parallel_solver")

#add_executable("pamg" pamg.cpp)
#target_link_dune_default_libraries("pamg")
add_subdirectory(pdelab)
#add_subdirectory(istl)
#add_subdirectory(super)
#add_subdirectory(heat_drp)
#add_subdirectory(heat_rp)
#add_subdirectory(setup_berlin)
#add_subdirectory(FE_vgl2)

add_subdirectory(FE_Newton_linear)
#add_subdirectory(FE_flame_propagation)
