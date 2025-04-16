LLVM_ENZYME=/usr/local/src/Enzyme/enzyme/build/Enzyme/LLVMEnzyme-16.so
FLAGS=-O3 -g -Wall -Wpedantic -Wconversion -Wsign-conversion -Wfloat-conversion

build: build/debug.exe build/descent.exe

run: run-debug run-descent

run-debug: build/debug.exe
	build/debug.exe
	convert real.ppm real.png
	convert gradient.ppm gradient.png

run-descent: build/descent.exe
	build/descent.exe
	convert -delay 10 -loop 1 $(shell echo temp/real_*.ppm) real.gif
	convert -delay 10 -loop 1 $(shell echo temp/gradient_*.ppm) gradient.gif

build/debug.exe: src/debug.cpp build/hello.o
	clang-16 $< $(FLAGS) -o $@

build/descent.exe: src/descent.cpp build/hello.o
	clang-16 $< $(FLAGS) -o $@

build/hello.o: build/output.ll
	clang-16 $< $(FLAGS) -lm -o $@

build/input.ll: src/hello.cpp
	clang-16 $< -S -emit-llvm -o $@ $(FLAGS)

build/output.ll: build/input.ll
	opt-16 $< --load-pass-plugin=$(LLVM_ENZYME) -passes=enzyme -o $@ -S

clean:
	rm -f build/*.ll build/*.exe build/*.o

.PHONY: build run run-descent run-debug clean
