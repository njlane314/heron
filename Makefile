CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra $(shell root-config --cflags)
LDFLAGS ?= $(shell root-config --libs) -lsqlite3

LIB_NAME = build/lib/libNuIO.so
LIB_SRC = lib/NuIO/src/ArtProvenanceIO.cxx \
          lib/NuIO/src/NuIOMaker.cxx \
          lib/NuIO/src/SampleIO.cxx
LIB_OBJ = $(LIB_SRC:.cxx=.o)

LIB_NUANA_NAME = build/lib/libNuAna.so
LIB_NUANA_SRC = lib/NuAna/src/SampleRDF.cxx
LIB_NUANA_OBJ = $(LIB_NUANA_SRC:.cxx=.o)

INCLUDES = -I./lib/NuIO/include -I./lib/NuAna/include

BIN_DIR = build/bin
NU_IO_MAKER = $(BIN_DIR)/nuIOmaker
BIN_SRC = bin/nuIOmaker/nuIOmaker.cxx
BIN_OBJ = $(BIN_SRC:.cxx=.o)

all: $(LIB_NAME) $(LIB_NUANA_NAME) $(NU_IO_MAKER)

$(LIB_NAME): $(LIB_OBJ)
	mkdir -p $(dir $(LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(LIB_OBJ) $(LDFLAGS) -o $(LIB_NAME)

$(LIB_NUANA_NAME): $(LIB_NUANA_OBJ)
	mkdir -p $(dir $(LIB_NUANA_NAME))
	$(CXX) -shared $(CXXFLAGS) $(LIB_NUANA_OBJ) $(LDFLAGS) -o $(LIB_NUANA_NAME)

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

$(NU_IO_MAKER): $(BIN_OBJ) $(LIB_OBJ)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(BIN_OBJ) $(LIB_OBJ) $(LDFLAGS) -o $(NU_IO_MAKER)

clean:
	rm -rf build
