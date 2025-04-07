# How to build

To build and run the demo, you should run

```sh
./make.sh run
```

This creates two images: `gradient.bpm` and `real.bpm`.

# References

This project uses [Enzyme](https://enzyme.mit.edu/) and [linmath.h](https://github.com/datenwolf/linmath.h).

Our dependency linmath.h is a simple header-only linear algebra library that defines types like `vec3`, `vec4`, and `mat4x4`.

Enzyme is a tool that autodiffs existing C code. The two most helpful pages in the reference are:

* <https://enzyme.mit.edu/getting_started/CallingConvention/>
* <https://enzyme.mit.edu/getting_started/Examples/>

Since we will be working with SDFs, Inigo Quilezles blog will be very helpful. There is a [page](https://iquilezles.org/articles/distfunctions/) that specifies SDFs for many 3d shapes, along with the operations needed to compose them to form more complex scenes.

I have also referenced Chris Wellons blog <https://nullprogram.com> in some places where we have used his subroutines.
