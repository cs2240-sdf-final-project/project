FLAGS=-Ofast -g -Wconversion -Wsign-conversion -Wfloat-conversion

build/a.exe: build/output.ll
	clang-16 $< $(FLAGS) -lm -o $@

build/input.ll: src/hello.c
	clang-16 $< -S -emit-llvm -o $@ $(FLAGS) 

build/output.ll: build/input.ll
	opt-16 $< --load-pass-plugin=/usr/local/src/Enzyme/enzyme/build/Enzyme/LLVMEnzyme-16.so -passes=enzyme -o $@ -S

clean:
	rm -f build/*.ll build/*.exe
