CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra $(shell root-config --cflags)
LDFLAGS ?= $(shell root-config --libs) -lsqlite3

LIB_NAME = build/lib/libNuIO.so
LIB_SRC = lib/NuIO/src/ArtProvenanceIO.cxx \
          lib/NuIO/src/NuIOMaker.cxx \
          lib/NuIO/src/SampleIO.cxx
LIB_OBJ = $(LIB_SRC:.cxx=.o)

INCLUDES = -I./lib/NuIO/include

BIN_DIR = build/bin
NU_IO_MAKER = $(BIN_DIR)/nuIOmaker
BIN_SRC = bin/nuIOmaker/nuIOmaker.cxx
BIN_OBJ = $(BIN_SRC:.cxx=.o)

all: $(LIB_NAME) $(NU_IO_MAKER)

$(LIB_NAME): $(LIB_OBJ)
	mkdir -p $(dir $(LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(LIB_OBJ) $(LDFLAGS) -o $(LIB_NAME)

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

$(NU_IO_MAKER): $(BIN_OBJ) $(LIB_OBJ)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(BIN_OBJ) $(LIB_OBJ) $(LDFLAGS) -o $(NU_IO_MAKER)

clean:
	rm -rf build
