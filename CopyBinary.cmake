if (binary_dir)
	if (EXE_NAME)
		set (TARGET_NAME ${EXE_NAME})
	endif (EXE_NAME)
	if (LIB_NAME)
		set (TARGET_NAME ${LIB_NAME})
	endif (LIB_NAME)
	
	ADD_CUSTOM_COMMAND(TARGET ${TARGET_NAME}
			  POST_BUILD
			  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${TARGET_NAME}> ${binary_dir}/.
	)
endif (binary_dir)

