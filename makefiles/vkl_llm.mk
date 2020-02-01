GPP = g++


############## START: VKL_LENS_LIGHT
CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  = -ljsoncpp -lvkl -lgfortran -lCCfits -lcfitsio -lgmp -lCGAL
EXT_LIB_DIR = common/vkl_lib/lib

ROOT_DIR  = lens_light_mass/vkl_llm
EXT_INC_DIR = common/vkl_lib/inc
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
	$(GPP) $(CPP_FLAGS) -I $(EXT_INC_DIR) -o $(BIN_DIR)/llm $(FULL_OBJ) $(CPP_LIBS) -L $(EXT_LIB_DIR) -Wl,-rpath,$(EXT_LIB_DIR)
lens_light_mass_clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*
############## END: VKL_LENS_LIGHT
