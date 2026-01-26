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

IO_LIB_NAME = build/lib/libNuxsecIO.so
IO_SRC = io/src/ArtFileProvenanceIO.cc \
         io/src/EventIO.cc \
         io/src/NormalisationService.cc \
         io/src/RunDatabaseService.cc \
         io/src/SampleIO.cc \
         io/src/SubRunInventoryService.cc
OBJ_DIR = build/obj
IO_OBJ = $(IO_SRC:%.cc=$(OBJ_DIR)/%.o)

ANA_LIB_NAME = build/lib/libNuxsecAna.so
ANA_SRC = ana/src/AnalysisConfigService.cc \
          ana/src/ColumnDerivationService.cc \
          ana/src/RDataFrameService.cc
ANA_OBJ = $(ANA_SRC:%.cc=$(OBJ_DIR)/%.o)

PLOT_LIB_NAME = build/lib/libNuxsecPlot.so
PLOT_SRC = plot/src/Plotter.cc \
           plot/src/StackedHist.cc
PLOT_OBJ = $(PLOT_SRC:%.cc=$(OBJ_DIR)/%.o)

NUXSEC_NAME = build/bin/nuxsec
NUXSEC_SRC = apps/src/nuxsec.cc

NUXSEC_ART_DRIVER_NAME = build/bin/nuxsecArtFileIOdriver
NUXSEC_ART_DRIVER_SRC = apps/src/nuxsecArtFileIOdriver.cc

NUXSEC_EVENT_DRIVER_NAME = build/bin/nuxsecEventIOdriver
NUXSEC_EVENT_DRIVER_SRC = apps/src/nuxsecEventIOdriver.cc

NUXSEC_SAMPLE_DRIVER_NAME = build/bin/nuxsecSampleIOdriver
NUXSEC_SAMPLE_DRIVER_SRC = apps/src/nuxsecSampleIOdriver.cc

INCLUDES = -I./io/include -I./ana/include -I./plot/include -I./apps/include

all: $(IO_LIB_NAME) $(ANA_LIB_NAME) $(PLOT_LIB_NAME) \
	 $(NUXSEC_NAME) $(NUXSEC_ART_DRIVER_NAME) $(NUXSEC_EVENT_DRIVER_NAME) \
	 $(NUXSEC_SAMPLE_DRIVER_NAME)

$(IO_LIB_NAME): $(IO_OBJ)
	mkdir -p $(dir $(IO_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(IO_OBJ) $(LDFLAGS) -o $(IO_LIB_NAME)

$(ANA_LIB_NAME): $(ANA_OBJ)
	mkdir -p $(dir $(ANA_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(ANA_OBJ) $(LDFLAGS) -o $(ANA_LIB_NAME)

$(PLOT_LIB_NAME): $(PLOT_OBJ)
	mkdir -p $(dir $(PLOT_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(PLOT_OBJ) $(LDFLAGS) -o $(PLOT_LIB_NAME)

$(NUXSEC_NAME): $(NUXSEC_SRC) $(IO_LIB_NAME) $(ANA_LIB_NAME)
	mkdir -p $(dir $(NUXSEC_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(NUXSEC_SRC) -Lbuild/lib -lNuxsecIO \
		-lNuxsecAna $(LDFLAGS) -o $(NUXSEC_NAME)

$(NUXSEC_ART_DRIVER_NAME): $(NUXSEC_ART_DRIVER_SRC) $(IO_LIB_NAME)
	mkdir -p $(dir $(NUXSEC_ART_DRIVER_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(NUXSEC_ART_DRIVER_SRC) -Lbuild/lib -lNuxsecIO \
		$(LDFLAGS) -o $(NUXSEC_ART_DRIVER_NAME)

$(NUXSEC_EVENT_DRIVER_NAME): $(NUXSEC_EVENT_DRIVER_SRC) $(IO_LIB_NAME) $(ANA_LIB_NAME)
	mkdir -p $(dir $(NUXSEC_EVENT_DRIVER_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(NUXSEC_EVENT_DRIVER_SRC) -Lbuild/lib -lNuxsecIO \
		-lNuxsecAna $(LDFLAGS) -o $(NUXSEC_EVENT_DRIVER_NAME)

$(NUXSEC_SAMPLE_DRIVER_NAME): $(NUXSEC_SAMPLE_DRIVER_SRC) $(IO_LIB_NAME)
	mkdir -p $(dir $(NUXSEC_SAMPLE_DRIVER_NAME))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(NUXSEC_SAMPLE_DRIVER_SRC) -Lbuild/lib -lNuxsecIO \
		$(LDFLAGS) -o $(NUXSEC_SAMPLE_DRIVER_NAME)

$(OBJ_DIR)/%.o: %.cc
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

clean:
	rm -rf build/lib build/bin build/obj

restart:
	rm -rf build/out
