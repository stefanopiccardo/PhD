add_executable(stabilization_test stabilization_test.cpp)
target_link_libraries(stabilization_test ${LINK_LIBS})
install(TARGETS stabilization_test RUNTIME DESTINATION bin)
