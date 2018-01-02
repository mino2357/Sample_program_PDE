CXX = g++
TARGET = main
CXXFLAGS = -std=c++1z -Wall -Wextra -O2 -march=native -mtune=native 
LDLFLAGS = -lstdc++
SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

run: all
	./$(TARGET)
clean:
	rm $(TARGET)
