# Add executables
add_executable(ShellOptimalDeformationNonLinAdaptive ShellOptimalDeformNonLinAdaptive.cpp)
add_executable(ShellOptimalDeformationIsometryLinAdaptive ShellOptimalDeformIsometryLinAdaptive.cpp)
add_executable(ShellOptimalDeformationIsometryNonLinAdaptive ShellOptimalDeformIsometryNonLinAdaptive.cpp)
# add_executable(ShellOptimalDeformationPlotter ShellOptimalDeformPlotter.cpp)
# add_executable(BashCompareThickness BashShellBucklingCompareThickness.cpp)

# Link the utility library
set(DKTELASTDEFOMRLIBS ${Boost_LIBRARIES} ${SUITESPARSE_LIBRARIES})

target_link_libraries(ShellOptimalDeformationNonLinAdaptive LINK_PUBLIC ${DKTELASTDEFOMRLIBS})
target_link_libraries(ShellOptimalDeformationIsometryLinAdaptive LINK_PUBLIC ${DKTELASTDEFOMRLIBS})
target_link_libraries(ShellOptimalDeformationIsometryNonLinAdaptive LINK_PUBLIC ${DKTELASTDEFOMRLIBS})
# target_link_libraries(ShellOptimalDeformationPlotter LINK_PUBLIC ${DKTELASTDEFOMRLIBS})
# target_link_libraries(BashCompareThickness LINK_PUBLIC ${DKTELASTDEFOMRLIBS})

if( PESOPT_WITH_IPOPT )
  target_link_libraries(ShellOptimalDeformationNonLinAdaptive LINK_PUBLIC ${IPOPT_LIBRARIES})
  target_link_libraries(ShellOptimalDeformationIsometryLinAdaptive LINK_PUBLIC ${IPOPT_LIBRARIES})
  target_link_libraries(ShellOptimalDeformationIsometryNonLinAdaptive LINK_PUBLIC ${IPOPT_LIBRARIES})
endif( PESOPT_WITH_IPOPT )

if( PESOPT_WITH_OPENMP )
  target_compile_definitions(ShellOptimalDeformationNonLinAdaptive PUBLIC PESOPT_WITH_OPENMP)
  set_target_properties(ShellOptimalDeformationNonLinAdaptive PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
  target_compile_definitions(ShellOptimalDeformationIsometryLinAdaptive PUBLIC PESOPT_WITH_OPENMP)
  set_target_properties(ShellOptimalDeformationIsometryLinAdaptive PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
  target_compile_definitions(ShellOptimalDeformationIsometryNonLinAdaptive PUBLIC PESOPT_WITH_OPENMP)
  set_target_properties(ShellOptimalDeformationIsometryNonLinAdaptive PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
#  target_compile_definitions(BashCompareThickness PUBLIC PESOPT_WITH_OPENMP)
#  set_target_properties(BashCompareThickness PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
endif( PESOPT_WITH_OPENMP )

if( PESOPT_WITH_VTK )
 target_link_libraries(ShellOptimalDeformationNonLinAdaptive LINK_PUBLIC ${VTK_LIBRARIES} )
 target_link_libraries(ShellOptimalDeformationIsometryLinAdaptive LINK_PUBLIC ${VTK_LIBRARIES} )
 target_link_libraries(ShellOptimalDeformationIsometryNonLinAdaptive LINK_PUBLIC ${VTK_LIBRARIES} )
#  target_link_libraries(ShellOptimalDeformationPlotter LINK_PUBLIC ${VTK_LIBRARIES} )
endif( PESOPT_WITH_VTK )
