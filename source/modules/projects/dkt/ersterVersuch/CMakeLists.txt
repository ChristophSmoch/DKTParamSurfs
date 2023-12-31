# Add executables

add_executable(ShellOptimalDeformationIsometrySemiNonLinAdaptive ShellOptimalDeformIsometrySemiNonLinAdaptive.cpp)
add_executable(TestEnergiesSemiNonLin TestEnergiesSemiNonLin.cpp)

# Link the utility library
set(DKTELASTDEFOMRSEMILIBS ${Boost_LIBRARIES} )


target_link_libraries(ShellOptimalDeformationIsometrySemiNonLinAdaptive LINK_PUBLIC ${DKTELASTDEFOMRSEMILIBS})
target_link_libraries(TestEnergiesSemiNonLin LINK_PUBLIC ${DKTELASTDEFOMRSEMILIBS})

if( PESOPT_WITH_IPOPT )
  target_link_libraries(ShellOptimalDeformationIsometrySemiNonLinAdaptive LINK_PUBLIC ${IPOPT_LIBRARIES})
  target_link_libraries(TestEnergiesSemiNonLin LINK_PUBLIC ${IPOPT_LIBRARIES})
endif( PESOPT_WITH_IPOPT )

if( PESOPT_WITH_OPENMP )
  target_compile_definitions(ShellOptimalDeformationIsometrySemiNonLinAdaptive PUBLIC PESOPT_WITH_OPENMP)
  set_target_properties(ShellOptimalDeformationIsometrySemiNonLinAdaptive PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
  target_compile_definitions(TestEnergiesSemiNonLin PUBLIC PESOPT_WITH_OPENMP)
  set_target_properties(TestEnergiesSemiNonLin PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
endif( PESOPT_WITH_OPENMP )

if( PESOPT_WITH_VTK )
  target_link_libraries(ShellOptimalDeformationIsometrySemiNonLinAdaptive LINK_PUBLIC ${VTK_LIBRARIES} )
  target_link_libraries(TestEnergiesSemiNonLin LINK_PUBLIC ${VTK_LIBRARIES} )
endif( PESOPT_WITH_VTK )
