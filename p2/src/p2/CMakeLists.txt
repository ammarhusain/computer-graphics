include_directories(
	${PROJECT_SOURCE_DIR}
	${SDL_INCLUDE_DIR}
)
add_executable(p2 main.cpp project.cpp)
target_link_libraries(p2 math application scene ${SDL_LIBRARY} ${PNG_LIBRARIES}
	${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${GLEW_LIBRARIES})

if(APPLE) 
	target_link_libraries(p2 SDLmain)
endif()

install(TARGETS p2 DESTINATION ${PROJECT_SOURCE_DIR}/..)
