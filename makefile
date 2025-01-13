CXX = /usr/bin/g++
TARGET = main
CXXFLAGS = -Wall -Wextra -O2 -march=native -mtune=native -std=c++2b 
LDLFLAGS = -lstdc++
SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

run: all
	./$(TARGET)

clean:
	rm $(TARGET)