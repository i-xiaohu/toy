cmake_minimum_required(VERSION 3.10)
project(mini_tools C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2 -Wno-unused-function")
add_definitions(-DUSE_MALLOC_WRAPPERS)

include_directories(.)

add_executable(mini_tools
		chr_extract.c
		hfastq.h
		hfastq.c
		kopen.c
		kseq.h
		ksort.h
		kstring.c
		kstring.h
		kthread.c
		kvec.h
		main.c
		malloc_wrap.c
		malloc_wrap.h
		pstat.c
		progress.c
		progress.h
		reads2fa.c
		samop.c
		samop.h
		table.c
		table.h
		utils.c
		utils.h
		wgsim_eval.c
		view_ref.c
		)

target_link_libraries(${PROJECT_NAME} -lz -lpthread)
if (CMAKE_HOST_UNIX)
	target_link_libraries(${PROJECT_NAME} -lrt)
endif ()