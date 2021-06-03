CXX = g++
INCLUDE = -I ../..
CXXFLAGS = -std=c++17 $(INCLUDE)
LDFLAGS = -lGLEW -lglfw -lGL 
SRCS = main.cpp ../../mathfunc/vec.cpp ../../mathfunc/matrix.cpp ../../opengl/visualize.cpp ../../opengl/window.cpp ../../opengl/vertarray.cpp ../../opengl/line.cpp  ../../opengl/point.cpp
OBJS = $(SRCS:.cpp=.o)
PROGRAM = Program

all: $(PROGRAM)
	./build/$(PROGRAM)

clean:
	rm -f $(OBJS)
	rm -f ./build/$(PROGRAM)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(PROGRAM) : $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o ./build/$(PROGRAM)
