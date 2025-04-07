LLVM_ENZYME=/usr/local/src/Enzyme/enzyme/build/Enzyme/LLVMEnzyme-16.so
FLAGS=-O3 -g -Wall -Wpedantic -Wconversion -Wsign-conversion -Wfloat-conversion

build: build/a.exe

run: build/a.exe
	build/a.exe

build/a.exe: build/output.ll
	clang-16 $< $(FLAGS) -lm -o $@

build/input.ll: src/hello.cpp
	clang-16 $< -S -emit-llvm -o $@ $(FLAGS) 

build/output.ll: build/input.ll
	opt-16 $< --load-pass-plugin=$(LLVM_ENZYME) -passes=enzyme -o $@ -S

clean:
	rm -f build/*.ll build/*.exe
