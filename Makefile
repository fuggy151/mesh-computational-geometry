CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -I./src/gKit

GKITSRC = $(wildcard src/gKit/*.cpp)

SRC = main.cpp $(GKITSRC)

OBJ = $(SRC:.cpp=.o)

TARGET = main

LDFLAGS = -lGL -lGLEW -lglfw -lSDL2 -lSDL2_image -ldl -lX11 -lpthread -lXrandr -lXi

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)