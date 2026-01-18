CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra $(shell root-config --cflags)
LDFLAGS ?= $(shell root-config --libs) -lsqlite3

IO_LIB_NAME = build/lib/libNuXsecIO.so
IO_SRC = io/src/ArtFileProvenanceRootIO.cc \
         io/src/SubrunTreeScanner.cc \
         io/src/RunInfoSqliteReader.cc \
         io/src/SampleTypes.cc \
         io/src/NuXSecTemplate.cc
IO_OBJ = $(IO_SRC:.cc=.o)

SAMPLE_LIB_NAME = build/lib/libNuXsecSample.so
SAMPLE_SRC = io/src/SampleAggregator.cc \
             io/src/SampleRootIO.cc
SAMPLE_OBJ = $(SAMPLE_SRC:.cc=.o)

ANA_LIB_NAME = build/lib/libNuXsecAna.so
ANA_SRC = ana/src/AnalysisRdfDefinitions.cc \
          ana/src/RDataFrameFactory.cc
ANA_OBJ = $(ANA_SRC:.cc=.o)

PLOT_LIB_NAME = build/lib/libNuXsecPlot.so
PLOT_SRC = plot/src/Plotter.cc \
           plot/src/StackedHist.cc
PLOT_OBJ = $(PLOT_SRC:.cc=.o)

RDF_BUILDER_NAME = build/bin/nuxsecSampleRDFbuilder
RDF_BUILDER_SRC = apps/src/nuxsecSampleRDFbuilder.cc

ART_AGGREGATOR_NAME = build/bin/nuxsecArtIOaggregator
ART_AGGREGATOR_SRC = apps/src/nuxsecArtIOaggregator.cc

SAMPLE_AGGREGATOR_NAME = build/bin/nuxsecSampleIOaggregator
SAMPLE_AGGREGATOR_SRC = apps/src/nuxsecSampleIOaggregator.cc

INCLUDES = -I./io/include -I./ana/include -I./plot/include -I./apps/include

all: $(IO_LIB_NAME) $(SAMPLE_LIB_NAME) $(ANA_LIB_NAME) $(PLOT_LIB_NAME) $(RDF_BUILDER_NAME) \
	 $(ART_AGGREGATOR_NAME) $(SAMPLE_AGGREGATOR_NAME)

$(IO_LIB_NAME): $(IO_OBJ)
	mkdir -p $(dir $(IO_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(IO_OBJ) $(LDFLAGS) -o $(IO_LIB_NAME)

$(SAMPLE_LIB_NAME): $(SAMPLE_OBJ) $(IO_LIB_NAME)
	mkdir -p $(dir $(SAMPLE_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(SAMPLE_OBJ) -Lbuild/lib -lNuXsecIO $(LDFLAGS) -o $(SAMPLE_LIB_NAME)

$(ANA_LIB_NAME): $(ANA_OBJ)
	mkdir -p $(dir $(ANA_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(ANA_OBJ) $(LDFLAGS) -o $(ANA_LIB_NAME)

$(PLOT_LIB_NAME): $(PLOT_OBJ)
	mkdir -p $(dir $(PLOT_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(PLOT_OBJ) $(LDFLAGS) -o $(PLOT_LIB_NAME)

$(RDF_BUILDER_NAME): $(RDF_BUILDER_SRC) $(IO_LIB_NAME) $(SAMPLE_LIB_NAME) $(ANA_LIB_NAME)
	mkdir -p $(dir $(RDF_BUILDER_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(RDF_BUILDER_SRC) -Lbuild/lib -lNuXsecIO -lNuXsecSample \
		-lNuXsecAna $(LDFLAGS) -o $(RDF_BUILDER_NAME)

$(ART_AGGREGATOR_NAME): $(ART_AGGREGATOR_SRC) $(IO_LIB_NAME)
	mkdir -p $(dir $(ART_AGGREGATOR_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(ART_AGGREGATOR_SRC) -Lbuild/lib -lNuXsecIO $(LDFLAGS) -o \
		$(ART_AGGREGATOR_NAME)

$(SAMPLE_AGGREGATOR_NAME): $(SAMPLE_AGGREGATOR_SRC) $(IO_LIB_NAME) $(SAMPLE_LIB_NAME)
	mkdir -p $(dir $(SAMPLE_AGGREGATOR_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SAMPLE_AGGREGATOR_SRC) -Lbuild/lib -lNuXsecSample -lNuXsecIO \
		$(LDFLAGS) -o $(SAMPLE_AGGREGATOR_NAME)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

clean:
	rm -rf build
