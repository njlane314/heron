CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra $(shell root-config --cflags)
LDFLAGS ?= $(shell root-config --libs) -lsqlite3

LIB_NAME = build/lib/libNuIO.so
LIB_SRC = src/lib/NuIO/src/StageResultIO.cxx
LIB_OBJ = $(LIB_SRC:.cxx=.o)

BIN_NAME = build/bin/nuIOaggregator.exe
BIN_SRC = src/bin/nuIOaggregator/nuIOaggregator.cxx
BIN_OBJ = $(BIN_SRC:.cxx=.o)

INCLUDES = -I./src/lib/NuIO/include

all: $(LIB_NAME) $(BIN_NAME)

$(LIB_NAME): $(LIB_OBJ)
	mkdir -p $(dir $(LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(LIB_OBJ) $(LDFLAGS) -o $(LIB_NAME)

$(BIN_NAME): $(BIN_OBJ) $(LIB_NAME)
	mkdir -p $(dir $(BIN_NAME))
	$(CXX) $(CXXFLAGS) $(BIN_OBJ) -Lbuild/lib -lNuIO $(LDFLAGS) -o $(BIN_NAME)

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

clean:
	rm -f $(LIB_OBJ) $(BIN_OBJ) $(LIB_NAME) $(BIN_NAME)
