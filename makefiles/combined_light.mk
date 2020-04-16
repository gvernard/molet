GPP = g++


CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  =  -lfftw3 -ljsoncpp -lvkl -lgfortran -lCCfits -lcfitsio -lgmp -lCGAL

ROOT_DIR = combined_light
SRC_DIR = $(ROOT_DIR)/src
INC_DIR = $(ROOT_DIR)/inc
BIN_DIR = $(ROOT_DIR)/bin
OBJ_DIR = $(ROOT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))


DEPS = mask_functions.hpp auxiliary_functions.hpp
OBJ  = mask_functions.o   auxiliary_functions.o   combine_light.o
FULL_DEPS = $(patsubst %,$(INC_DIR)/%,$(DEPS)) #Pad names with dir
FULL_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir
#$(info $$OBJ is [${FULL_DEPS}])

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(FULL_DEPS)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -c -o $@ $<

combined_light: $(FULL_OBJ)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -o $(BIN_DIR)/combine_light $(FULL_OBJ) $(CPP_LIBS)
combined_light_clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*

