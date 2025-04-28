LLVM_ENZYME=/usr/local/src/Enzyme/enzyme/build/Enzyme/LLVMEnzyme-16.so
FLAGS=-O2 -g -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wfloat-conversion -fopenmp
FRONTEND_FLAGS=-fsanitize=address,undefined
LINK_FLAGS=-lm -lpthread

build: build/debug.exe build/descent.exe

run: run-debug run-descent

run-debug: build/debug.exe
	mkdir -p debug-gradient
	rm -f debug-gradient/*
	build/debug.exe

run-descent: build/descent.exe
	mkdir -p descent-sequence
	rm -f descent-sequence/*
	build/descent.exe
	ffmpeg -y -framerate 30 -pattern_type glob -i 'descent-sequence/real_*.ppm' descent-sequence.mp4

build/debug.exe: src/debug.cpp build/hello.o
	clang++-16 $^ $(FRONTEND_FLAGS) $(FLAGS) $(LINK_FLAGS) -o $@

build/descent.exe: src/descent.cpp build/hello.o
	clang++-16 $^ $(FRONTEND_FLAGS) $(FLAGS) $(LINK_FLAGS) -o $@

build/hello.o: build/output.ll
	clang++-16 -c $^ $(FRONTEND_FLAGS) $(FLAGS) -o $@

build/input.ll: src/hello.cpp
	clang++-16 $^ -S -emit-llvm -o $@ $(FLAGS)

build/output.ll: build/input.ll
	opt-16 $^ --load-pass-plugin=$(LLVM_ENZYME) -passes=enzyme -o $@ -S

clean:
	rm -f build/*.ll build/*.exe build/*.o

.PHONY: build run run-descent run-debug clean
