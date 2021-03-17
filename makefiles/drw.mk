.DEFAULT_GOAL := drw

GPP = g++


CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  = -lfftw3 -ljsoncpp -lgerlumph -lpng -lCCfits -lcfitsio
EXT_LIBS  = -linstruments

EXT_LIB_DIR = instrument_modules/lib
EXT_INC_DIR = instrument_modules/include

ROOT_DIR = variability/intrinsic/DRW
SRC_DIR = $(ROOT_DIR)/src
INC_DIR = $(ROOT_DIR)/inc
BIN_DIR = $(ROOT_DIR)/bin
OBJ_DIR = $(ROOT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))


HEADERS = $(shell find $(INC_DIR) -type f -name '*.hpp')
OBJ  = auxiliary_functions.o drw.o
FULL_DEPS = $(patsubst %,$(INC_DIR)/%,$(DEPS)) #Pad names with dir
FULL_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir

#$(info $$OBJ is [${HEADERS}])
#$(info $$OBJ is [${FULL_OBJ}])


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(FULL_DEPS)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -I $(EXT_INC_DIR) -c -o $@ $<

drw: $(FULL_OBJ)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -I $(EXT_INC_DIR) -o $(BIN_DIR)/drw $(FULL_OBJ) $(CPP_LIBS) $(EXT_LIBS) -L $(EXT_LIB_DIR) -Wl,-rpath,$(EXT_LIB_DIR)
clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*

