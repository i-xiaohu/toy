cmake_minimum_required(VERSION 3.13)
project(toy C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(toy
		check2Files.c
		chr_extract.c
		eva2sam.c
		generate_cs.c
		hfastq.c
		hsfq.c
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
		parallel_test.c
		proc_stat.c
		progress.c
		progress.h
		reads2fa.c
		sam2sfq.c
		samop.c
		samop.h
		sync_pe.c
		sync_pe.h
		table.c
		table.h
		utils.c
		utils.h)
