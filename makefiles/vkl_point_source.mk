.DEFAULT_GOAL := point_source

GPP = g++


CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  =  -ljsoncpp -lvkl -lgfortran -lCCfits -lcfitsio -lgmp -lCGAL
EXT_LIBS  = -linstruments

EXT_LIB_DIR = instrument_modules/lib
EXT_INC_DIR = instrument_modules/include

ROOT_DIR = lensed_point_source/vkl_point_source
SRC_DIR = $(ROOT_DIR)/src
INC_DIR = $(ROOT_DIR)/inc
BIN_DIR = $(ROOT_DIR)/bin
OBJ_DIR = $(ROOT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))


HEADERS = $(shell find $(INC_DIR) -type f -name '*.hpp')
OBJ  = polygons.o   point_source.o
FULL_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir
#$(info $$OBJ is [${FULL_DEPS}])

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -I $(EXT_INC_DIR) -c -o $@ $<

point_source: $(FULL_OBJ)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -I $(EXT_INC_DIR) -o $(BIN_DIR)/point_source $(FULL_OBJ) $(CPP_LIBS) $(EXT_LIBS) -L $(EXT_LIB_DIR) -Wl,-rpath,$(EXT_LIB_DIR)
clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*
