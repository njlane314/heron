CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra $(shell root-config --cflags)
LDFLAGS ?= $(shell root-config --libs) -lsqlite3

LIB_NAME = libNuIO.so
LIB_SRC = src/lib/NuIO/src/StageResultIO.cxx
LIB_OBJ = $(LIB_SRC:.cxx=.o)

BIN_NAME = nuIOcondenser.exe
BIN_SRC = src/bin/nuIOcondenser/nuIOcondenser.cxx
BIN_OBJ = $(BIN_SRC:.cxx=.o)

INCLUDES = -I./src/lib/NuIO/include

all: $(LIB_NAME) $(BIN_NAME)

$(LIB_NAME): $(LIB_OBJ)
	$(CXX) -shared $(CXXFLAGS) $(LIB_OBJ) $(LDFLAGS) -o $(LIB_NAME)

$(BIN_NAME): $(BIN_OBJ) $(LIB_NAME)
	$(CXX) $(CXXFLAGS) $(BIN_OBJ) -L. -lNuIO $(LDFLAGS) -o $(BIN_NAME)

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

clean:
	rm -f $(LIB_OBJ) $(BIN_OBJ) $(LIB_NAME) $(BIN_NAME)
