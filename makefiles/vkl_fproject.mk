GPP = g++


############## START: VKL_FPROJECT
CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  = -ljsoncpp -lvkl -lgfortran -lCCfits -lcfitsio -lgmp -lCGAL
EXT_LIB_DIR = common/vkl_lib/lib

ROOT_DIR = lensed_extended_source/vkl_fproject
EXT_INC_DIR = common/vkl_lib/inc
SRC_DIR = $(ROOT_DIR)/src
BIN_DIR = $(ROOT_DIR)/bin
OBJ_DIR = $(ROOT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))


OBJ  = fproject.o
FULL_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir
#$(info $$OBJ is [${FULL_DEPS}])

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(GPP) $(CPP_FLAGS) -I $(EXT_INC_DIR) -c -o $@ $<

fproject: $(FULL_OBJ)
	$(GPP) $(CPP_FLAGS) -I $(EXT_INC_DIR) -o $(BIN_DIR)/fproject $(FULL_OBJ) $(CPP_LIBS) -L $(EXT_LIB_DIR) -Wl,-rpath,$(EXT_LIB_DIR)
fproject_clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*
############## END: VKL_FPROJECT
