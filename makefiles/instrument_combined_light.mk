GPP = g++


############## START: Fproject
CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  =  -ljsoncpp -lvkl -lgfortran -lCCfits -lcfitsio -lgmp -lCGAL
EXT_LIB_DIR = common/vkl_lib/lib

ROOT_DIR = instrument_combined_light
EXT_INC_DIR = common/vkl_lib/inc
SRC_DIR = $(ROOT_DIR)/src
INC_DIR = $(ROOT_DIR)/inc
BIN_DIR = $(ROOT_DIR)/bin
OBJ_DIR = $(ROOT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))


DEPS = auxiliary_functions.hpp
OBJ  = auxiliary_functions.o   combine_light.o
FULL_DEPS = $(patsubst %,$(INC_DIR)/%,$(DEPS)) #Pad names with dir
FULL_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir
#$(info $$OBJ is [${FULL_DEPS}])

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(FULL_DEPS)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -I $(EXT_INC_DIR) -c -o $@ $<

instrument_combined_light: $(FULL_OBJ)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -I $(EXT_INC_DIR) -o $(BIN_DIR)/combine_light $(FULL_OBJ) $(CPP_LIBS) -L $(EXT_LIB_DIR) -Wl,-rpath,$(EXT_LIB_DIR)
instrument_combined_light_clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*
############## END: Fproject
