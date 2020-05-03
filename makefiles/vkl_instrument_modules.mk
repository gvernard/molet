# The code MUST be compilled by providing INSTRUMENT_PATH=/path/to/instrument/files/without/quotes/and/ending/with/slash/ to make
# To specify a different path you have to call "make clean" first.
.DEFAULT_GOAL := instrument_modules

GPP = g++

CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  = -lvkl -lfftw3 -ljsoncpp

ROOT_DIR = instrument_modules
SRC_DIR = $(ROOT_DIR)/src
INC_DIR = $(ROOT_DIR)/include
LIB_DIR = $(ROOT_DIR)/lib
OBJ_DIR = $(ROOT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(LIB_DIR))


HEADERS = $(shell find $(HEADER_DIR) -type f -name '*.hpp')

SPECIAL = instruments.cpp
FULL_SPECIAL = $(patsubst %, $(SRC_DIR)/%,$(SPECIAL))
OBJ_SPECIAL  = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(FULL_SPECIAL:.cpp=.o))

SOURCES = noise.cpp
FULL_SOURCES = $(patsubst %, $(SRC_DIR)/%,$(SOURCES))
OBJ_SOURCES  = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(FULL_SOURCES:.cpp=.o))


#$(info $$OBJ is [${FULL_SPECIAL}])
#$(info $$OBJ is [${OBJ_SPECIAL}])


QUOTED_INSTRUMENT_PATH = $(addprefix '",$(addsuffix "',$(INSTRUMENT_PATH)))
INSTRUMENT_PATH_FLAGS = -DINSTRUMENT_PATH=$(QUOTED_INSTRUMENT_PATH)

# INSTRUMENT PATH
$(OBJ_DIR)/instruments.o: $(FULL_SPECIAL) $(HEADERS)
ifndef INSTRUMENT_PATH
	$(error INSTRUMENT_PATH is not set)
else
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) $(INSTRUMENT_PATH_FLAGS) -c -o $@ $<
endif

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -c -o $@ $<

instrument_modules: $(OBJ_SPECIAL) $(OBJ_SOURCES)
	$(GPP) -shared -Wl,-soname,libinstruments.so -o $(LIB_DIR)/libinstruments.so $(OBJ_SPECIAL) $(OBJ_SOURCES) $(CPP_LIBS)
clean:
	$(RM) -r $(OBJ_DIR)/* $(LIB_DIR)/*
