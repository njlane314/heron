CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra $(shell root-config --cflags)
LDFLAGS ?= $(shell root-config --libs) -lsqlite3

LIB_NAME = build/lib/libNuIO.so
LIB_SRC = lib/NuIO/src/ArtProvenanceIO.cxx \
          lib/NuIO/src/SampleIO.cxx
LIB_OBJ = $(LIB_SRC:.cxx=.o)

INCLUDES = -I./lib/NuIO/include

all: $(LIB_NAME)

$(LIB_NAME): $(LIB_OBJ)
	mkdir -p $(dir $(LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(LIB_OBJ) $(LDFLAGS) -o $(LIB_NAME)

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

clean:
	rm -rf build
