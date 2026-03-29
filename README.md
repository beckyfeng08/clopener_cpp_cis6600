# Clopening
## Running it
After any edits, rebuild the project and run it from the `Build` dir. Without any command line args it defaults to the bunny mesh but you can use another trimesh by adding a local filepath as a command line arg. (Add your own trimeshes to the `data` directory.) 

``` .\Release\example.exe ../data/eiffel.obj```

Once the project is running, press the `c` key to perform closing.
## Parameters
`closing_flow.h` has the params for the closing flow function. Note that I've been using `max_iter=1`; change this for testing the full clopening loop once we have the curvature flow step working. I think we will also have to do some tuning for the `h` parameter (target edge length during remeshing). Mayhaps it should be adaptive.

## Project to-do list
- Enable non triangle meshes (currently assume trimesh)
- Curvature flow update over time step


=============================================================================
# libigl example project

A blank project example showing how to use libigl and cmake. Feel free and
encouraged to copy or fork this project as a way of starting a new personal
project using libigl.

## Compile

Compile this project using the standard cmake routine:

```
mkdir build
cd build
cmake ..
make
```

This should find and build the dependencies and create a `example` binary.

For Windows, use 

`cmake --build . --config Release`  
instead of the last line `make`. This will create a directory `Release` containing an executable,`example.exe`.

## Run

From within the `build` directory just issue:

```
./example
```

A glfw app should launch displaying a 3D cube.

For Windows, run

`.\Release\example.exe`.

## Using other modules of libigl

This example project uses the `igl::opengl::glfw::Viewer`, therefore it requires
the glfw module of libigl. This shows up in the CMakeLists.txt 

```cmake
igl_include(glfw)
…
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)
```

Suppose you also wanted to use the triangle module in libigl. Then you would
change these to

```cmake
igl_include(glfw)
igl_include(restricted triangle)
…
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw igl_restricted::triangle)
```

The "restricted" appears in this case because the triangle library has a more
restricted license than libigl. See other examples commented out in
CMakeLists.txt.

## Dependencies

The only dependencies are STL, Eigen, [libigl](http://libigl.github.io/libigl/) and the dependencies
of the `igl::opengl::glfw::Viewer` (OpenGL, glad and GLFW).

The CMake build system will automatically download libigl and its dependencies using
[CMake FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html),
thus requiring no setup on your part.

### Use a local copy of libigl

You can use the CMake cache variable `FETCHCONTENT_SOURCE_DIR_LIBIGL` when configuring your CMake project for
the first time to aim it at a local copy of libigl instead.

```
cmake -DFETCHCONTENT_SOURCE_DIR_LIBIGL=<path-to-libigl> ..
```

When changing this value, do not forget to clear your `CMakeCache.txt`, or to update the cache variable
via `cmake-gui` or `ccmake`.
