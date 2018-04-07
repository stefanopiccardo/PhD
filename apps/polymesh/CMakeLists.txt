add_executable(polymesh polymesh.cpp)
target_link_libraries(polymesh ${LINK_LIBS})
install(TARGETS polymesh RUNTIME DESTINATION bin)
