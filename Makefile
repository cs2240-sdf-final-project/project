build/a.exe: build/output.ll
	clang-16 $< -O3 -o $@

build/input.ll: src/hello.c
	clang-16 $< -S -I3rd/linmath -emit-llvm -o $@ -Ofast -fno-vectorize -fno-slp-vectorize -fno-unroll-loops

build/output.ll: build/input.ll
	opt-16 $< --load-pass-plugin=/usr/local/src/Enzyme/enzyme/build/Enzyme/LLVMEnzyme-16.so -passes=enzyme -o $@ -S

clean:
	rm -f build/*.ll build/*.exe
