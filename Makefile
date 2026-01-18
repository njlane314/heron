CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra $(shell root-config --cflags)
LDFLAGS ?= $(shell root-config --libs) -lsqlite3

IO_LIB_NAME = build/lib/libNuxsecIO.so
IO_SRC = io/src/ArtProvenanceIO.cc \
         io/src/SampleIO.cc \
         io/src/SubRunScanner.cc \
         io/src/RunInfoDB.cc
IO_OBJ = $(IO_SRC:.cc=.o)

ANA_LIB_NAME = build/lib/libNuxsecAna.so
ANA_SRC = ana/src/AnalysisProcessor.cc \
          ana/src/RDFBuilder.cc
ANA_OBJ = $(ANA_SRC:.cc=.o)

RDF_BUILDER_NAME = build/bin/sampleRDFbuilder
RDF_BUILDER_SRC = apps/src/sampleRDFbuilder.cc

ART_AGGREGATOR_NAME = build/bin/artIOaggregator
ART_AGGREGATOR_SRC = apps/src/artIOaggregator.cc

SAMPLE_AGGREGATOR_NAME = build/bin/sampleIOaggregator
SAMPLE_AGGREGATOR_SRC = apps/src/sampleIOaggregator.cc

INCLUDES = -I./io/include -I./ana/include -I./apps/include

all: $(IO_LIB_NAME) $(ANA_LIB_NAME) $(RDF_BUILDER_NAME) $(ART_AGGREGATOR_NAME) \
	 $(SAMPLE_AGGREGATOR_NAME)

$(IO_LIB_NAME): $(IO_OBJ)
	mkdir -p $(dir $(IO_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(IO_OBJ) $(LDFLAGS) -o $(IO_LIB_NAME)

$(ANA_LIB_NAME): $(ANA_OBJ)
	mkdir -p $(dir $(ANA_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(ANA_OBJ) $(LDFLAGS) -o $(ANA_LIB_NAME)

$(RDF_BUILDER_NAME): $(RDF_BUILDER_SRC) $(IO_LIB_NAME) $(ANA_LIB_NAME)
	mkdir -p $(dir $(RDF_BUILDER_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(RDF_BUILDER_SRC) -Lbuild/lib -lNuxsecIO -lNuxsecAna \
		$(LDFLAGS) -o $(RDF_BUILDER_NAME)

$(ART_AGGREGATOR_NAME): $(ART_AGGREGATOR_SRC) $(IO_LIB_NAME)
	mkdir -p $(dir $(ART_AGGREGATOR_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(ART_AGGREGATOR_SRC) -Lbuild/lib -lNuxsecIO $(LDFLAGS) -o \
		$(ART_AGGREGATOR_NAME)

$(SAMPLE_AGGREGATOR_NAME): $(SAMPLE_AGGREGATOR_SRC) $(IO_LIB_NAME)
	mkdir -p $(dir $(SAMPLE_AGGREGATOR_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SAMPLE_AGGREGATOR_SRC) -Lbuild/lib -lNuxsecIO $(LDFLAGS) -o \
		$(SAMPLE_AGGREGATOR_NAME)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

clean:
	rm -rf build
