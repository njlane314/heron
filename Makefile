CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra $(shell root-config --cflags)
LDFLAGS ?= $(shell root-config --libs) -lsqlite3

LIB_NAME = build/lib/libNuIO.so
LIB_SRC = lib/NuIO/src/ArtProvenanceIO.cxx \
          lib/NuIO/src/SampleIO.cxx
LIB_OBJ = $(LIB_SRC:.cxx=.o)

LIB_NUANA_NAME = build/lib/libNuAna.so
LIB_NUANA_SRC = lib/NuAna/src/SampleRDF.cxx
LIB_NUANA_OBJ = $(LIB_NUANA_SRC:.cxx=.o)

RDF_BUILDER_NAME = build/bin/sampleRDFbuilder
RDF_BUILDER_SRC = bin/sampleRDFbuilder/sampleRDFbuilder.cxx

INCLUDES = -I./lib/NuIO/include -I./lib/NuAna/include

all: $(LIB_NAME) $(LIB_NUANA_NAME) $(RDF_BUILDER_NAME)

$(LIB_NAME): $(LIB_OBJ)
	mkdir -p $(dir $(LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(LIB_OBJ) $(LDFLAGS) -o $(LIB_NAME)

$(LIB_NUANA_NAME): $(LIB_NUANA_OBJ)
	mkdir -p $(dir $(LIB_NUANA_NAME))
	$(CXX) -shared $(CXXFLAGS) $(LIB_NUANA_OBJ) $(LDFLAGS) -o $(LIB_NUANA_NAME)

$(RDF_BUILDER_NAME): $(RDF_BUILDER_SRC) $(LIB_NAME) $(LIB_NUANA_NAME)
	mkdir -p $(dir $(RDF_BUILDER_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(RDF_BUILDER_SRC) -Lbuild/lib -lNuIO -lNuAna $(LDFLAGS) -o $(RDF_BUILDER_NAME)

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

clean:
	rm -rf build
