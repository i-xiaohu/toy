cmake_minimum_required(VERSION 3.13)
project(toy C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(toy
		chr_extract.c
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
		proc_stat.c
		progress.c
		progress.h
		reads2fa.c
		sam2sfq.c
		samop.c
		samop.h
		table.c
		table.h
		utils.c
		utils.h
		wgsim_eval.c
		view_ref.c hfastq.h)
