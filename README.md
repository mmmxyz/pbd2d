# pbd2d

requirements
- GLFW
- GLEW

例:project/freefallを実行する
```
$ cd pbd2d/project/freefall
$ mkdir build
$ g++ -std=c++17 -I ../.. -c main.cpp -o main.o
$ g++ -std=c++17 -I ../.. -c ../../mathfunc/vec.cpp -o ../../mathfunc/vec.o
$ g++ -std=c++17 -I ../.. -c ../../mathfunc/matrix.cpp -o ../../mathfunc/matrix.o
$ g++ -std=c++17 -I ../.. -c ../../opengl/visualize.cpp -o ../../opengl/visualize.o
$ g++ -std=c++17 -I ../.. -c ../../opengl/window.cpp -o ../../opengl/window.o
$ g++ -std=c++17 -I ../.. -c ../../opengl/vertarray.cpp -o ../../opengl/vertarray.o
$ g++ -std=c++17 -I ../.. -c ../../opengl/line.cpp -o ../../opengl/line.o
$ g++ -std=c++17 -I ../.. -c ../../opengl/point.cpp -o ../../opengl/point.o
$ g++ main.o ../../mathfunc/vec.o ../../mathfunc/matrix.o ../../opengl/visualize.o ../../opengl/window.o ../../opengl/vertarray.o ../../opengl/line.o ../../opengl/point.o -lGLEW -lglfw -lGL  -o ./build/Program
$ ./build/Program
```
あるいはCMakeで
```
$ cd pbd2d/project/freefall
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build .
````
