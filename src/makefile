CC = gcc
CXX = g++
DEBUG = -g
CXXFLAGS = $(shell fltk-config --use-gl --use-images --cxxflags ) -I.
LDFLAGS = $(shell fltk-config --use-gl --use-glut --use-images --ldflags )
LDSTATIC = $(shell fltk-config --use-gl --use-glut --use-images --ldstaticflags )

TARGET = runMaze buildMaze
SRCS = Cell.cpp Edge.cpp LineSeg.cpp Vector3.cpp \
    MapWidget.cpp Maze.cpp MazeWindow.cpp RunMaze.cpp Vertex.cpp
OBJS = $(SRCS:.cpp=.o)

BUILDSRCS = Cell.cpp Edge.cpp LineSeg.cpp Vector3.cpp \
	MapWidget.cpp Maze.cpp MazeWindow.cpp BuildMaze.cpp Vertex.cpp
BUILDOBJS=$(BUILDSRCS:.cpp=.o)

.SUFFIXES: .o .cpp
*.o: *.cpp
	$(CXX) $(CXXFLAGS) $(DEBUG) -c $<

all: $(TARGET)
	# $(CXX) -o $(TARGET) $(OBJS) $(LDSTATIC)

runmaze: $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS)

buildmaze: $(BUILDOBJS)
	$(CXX) -o $@ $(BUILDOBJS) $(LDFLAGS)

depend:
	makedepend -- $(INCLUDE) -- $(SRCS) $(BUILDSRCS)

clean:
	rm -f *.o 2> /dev/null
	rm -f $(TARGET) 2> /dev/null

# $(TARGET): $(OBJS)
