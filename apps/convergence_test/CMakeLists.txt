add_executable(convergence_test convergence_test.cpp)
target_link_libraries(convergence_test ${LINK_LIBS})
install(TARGETS convergence_test RUNTIME DESTINATION bin)
