CXX ?= $(shell $(ROOT_CONFIG) --cxx)
ROOT_CONFIG ?= root-config

ifeq ($(shell command -v $(ROOT_CONFIG) 2>/dev/null),)
$(error ROOT not found. Please set up ROOT so that 'root-config' is on your PATH.)
endif

NLOHMANN_JSON_INC ?=
NLOHMANN_JSON_CFLAGS ?=
ifneq ($(strip $(NLOHMANN_JSON_INC)),)
NLOHMANN_JSON_CFLAGS := -isystem $(NLOHMANN_JSON_INC)
endif

CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra $(shell $(ROOT_CONFIG) --cflags) $(NLOHMANN_JSON_CFLAGS)
LDFLAGS ?= $(shell $(ROOT_CONFIG) --libs) -lsqlite3

IO_LIB_NAME = build/lib/libNuXsecIO.so
IO_SRC = io/src/ArtFileProvenanceRootIO.cc \
         io/src/SubrunTreeScanner.cc \
         io/src/RunInfoSqliteReader.cc \
         io/src/SampleTypes.cc \
         io/src/NuXSecTemplate.cc \
         io/src/TemplateRootIO.cc
IO_OBJ = $(IO_SRC:.cc=.o)

SAMPLE_LIB_NAME = build/lib/libNuXsecSample.so
SAMPLE_SRC = io/src/SampleAggregator.cc \
             io/src/SampleRootIO.cc
SAMPLE_OBJ = $(SAMPLE_SRC:.cc=.o)

ANA_LIB_NAME = build/lib/libNuXsecAna.so
ANA_SRC = ana/src/AnalysisRdfDefinitions.cc \
          ana/src/RDataFrameFactory.cc \
          ana/src/TemplateSpec.cc
ANA_OBJ = $(ANA_SRC:.cc=.o)

PLOT_LIB_NAME = build/lib/libNuXsecPlot.so
PLOT_SRC = plot/src/Plotter.cc \
           plot/src/StackedHist.cc
PLOT_OBJ = $(PLOT_SRC:.cc=.o)

ART_AGGREGATOR_NAME = build/bin/nuxsecArtIOaggregator
ART_AGGREGATOR_SRC = apps/src/nuxsecArtIOaggregator.cc

SAMPLE_AGGREGATOR_NAME = build/bin/nuxsecSampleIOaggregator
SAMPLE_AGGREGATOR_SRC = apps/src/nuxsecSampleIOaggregator.cc

TEMPLATE_MAKER_NAME = build/bin/nuxsecTemplateMaker
TEMPLATE_MAKER_SRC = apps/src/nuxsecTemplateMaker.cc

INCLUDES = -I./io/include -I./ana/include -I./plot/include -I./apps/include

all: $(IO_LIB_NAME) $(SAMPLE_LIB_NAME) $(ANA_LIB_NAME) $(PLOT_LIB_NAME) \
	 $(ART_AGGREGATOR_NAME) $(SAMPLE_AGGREGATOR_NAME) $(TEMPLATE_MAKER_NAME)

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

$(ART_AGGREGATOR_NAME): $(ART_AGGREGATOR_SRC) $(IO_LIB_NAME)
	mkdir -p $(dir $(ART_AGGREGATOR_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(ART_AGGREGATOR_SRC) -Lbuild/lib -lNuXsecIO $(LDFLAGS) -o \
		$(ART_AGGREGATOR_NAME)

$(SAMPLE_AGGREGATOR_NAME): $(SAMPLE_AGGREGATOR_SRC) $(IO_LIB_NAME) $(SAMPLE_LIB_NAME)
	mkdir -p $(dir $(SAMPLE_AGGREGATOR_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SAMPLE_AGGREGATOR_SRC) -Lbuild/lib -lNuXsecSample -lNuXsecIO \
		$(LDFLAGS) -o $(SAMPLE_AGGREGATOR_NAME)

$(TEMPLATE_MAKER_NAME): $(TEMPLATE_MAKER_SRC) $(IO_LIB_NAME) $(SAMPLE_LIB_NAME) $(ANA_LIB_NAME)
	mkdir -p $(dir $(TEMPLATE_MAKER_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(TEMPLATE_MAKER_SRC) -Lbuild/lib -lNuXsecSample -lNuXsecIO \
		-lNuXsecAna $(LDFLAGS) -o $(TEMPLATE_MAKER_NAME)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

clean:
	rm -rf build
