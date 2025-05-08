LLVM_ENZYME=/usr/local/src/Enzyme/enzyme/build/Enzyme/LLVMEnzyme-16.so
FLAGS=-O3 -g -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wfloat-conversion -fopenmp
FRONTEND_FLAGS=-fsanitize=address,undefined
LINK_FLAGS=-lm -lpthread

build: build/debug.exe build/descent.exe

run: run-debug run-descent

run-debug: build/debug.exe
	build/debug.exe

run-descent: build/descent.exe
	rm -f descent-sequence/*.ppm
	build/descent.exe
	ffmpeg -y -framerate 30 -pattern_type glob -i 'descent-sequence/real_*.ppm' descent-sequence.mp4

build/debug.exe: src/debug.cpp build/render.o build/sim_random.o build/image.o
	clang++-16 $^ $(FRONTEND_FLAGS) $(FLAGS) $(LINK_FLAGS) -o $@

build/descent.exe: src/descent.cpp build/render.o build/sim_random.o build/image.o
	clang++-16 $^ $(FRONTEND_FLAGS) $(FLAGS) $(LINK_FLAGS) -o $@

build/render.o: build/render.2.ll
	clang++-16 -c $^ $(FRONTEND_FLAGS) $(FLAGS) -o $@

build/render.2.ll: build/render.1.ll
	opt-16 $^ --load-pass-plugin=$(LLVM_ENZYME) '-passes=enzyme' -o $@ -S

build/render.1.ll: src/render.cpp
	clang++-16 $^ -S -emit-llvm -o $@ $(FLAGS)

build/sim_random.o: src/sim_random.cpp
	clang++-16 -c $^ $(FLAGS) $(FRONTEND_FLAGS) -o $@

build/image.o: src/image.cpp
	clang++-16 -c $^ $(FLAGS) $(FRONTEND_FLAGS) -o $@

clean:
	rm -f build/*.ll build/*.exe build/*.o

.PHONY: build run run-descent run-debug clean
