#add_executable(test1 test1.cpp)
#target_link_libraries(test1 ${LINK_LIBS})


find_package(Lua)
if (LUA_FOUND)
    include_directories("${LUA_INCLUDE_DIR}")
    set(LINK_LIBS ${LINK_LIBS} ${LUA_LIBRARIES})

    add_executable(cuthho_square cuthho_square.cpp)
    target_link_libraries(cuthho_square ${LINK_LIBS})
endif ()

add_executable(convergence_test convergence_test.cpp)
target_link_libraries(convergence_test ${LINK_LIBS})
