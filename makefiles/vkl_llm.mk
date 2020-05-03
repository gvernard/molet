.DEFAULT_GOAL := lens_light_mass

GPP = g++



CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  = -ljsoncpp -lvkl -lgfortran -lCCfits -lcfitsio -lgmp -lCGAL
EXT_LIBS  = -linstruments

EXT_LIB_DIR = instrument_modules/lib
EXT_INC_DIR = instrument_modules/include

ROOT_DIR  = lens_light_mass/vkl_llm
SRC_DIR = $(ROOT_DIR)/src
BIN_DIR = $(ROOT_DIR)/bin
OBJ_DIR = $(ROOT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))

OBJ  = lens_light_mass.o
FULL_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir
#$(info $$OBJ is [${LL_OBJ}])


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(GPP) $(CPP_FLAGS) -I $(EXT_INC_DIR) -c -o $@ $<

lens_light_mass: $(FULL_OBJ)
	$(GPP) $(CPP_FLAGS) -I $(EXT_INC_DIR) -o $(BIN_DIR)/llm $(FULL_OBJ) $(CPP_LIBS) $(EXT_LIBS) -L $(EXT_LIB_DIR) -Wl,-rpath,$(EXT_LIB_DIR)
clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*
