GPP = g++


############## START: vkl_lens_light
LL_FLAGS = -std=c++11 -fPIC -g -frounding-math
LL_LIBS  = -lCCfits -lcfitsio -ljsoncpp

LL_DIR  = lens_light/vkl_lens_light
SRC_DIR = $(LL_DIR)/src
INC_DIR = $(LL_DIR)/inc
BIN_DIR = $(LL_DIR)/bin
OBJ_DIR = $(LL_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))

DEPS = imagePlane.hpp sourceProfile.hpp tableDefinition.hpp
OBJ  = imagePlane.o   sourceProfile.o   lens_light.o
LL_DEPS = $(patsubst %,$(INC_DIR)/%,$(DEPS)) #Pad names with dir
LL_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir
#$(info $$OBJ is [${LL_OBJ}])


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(LL_DEPS)
	$(GPP) $(LL_FLAGS) -I $(INC_DIR) -c -o $@ $<
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(FORTRAN) $(FORTRAN_FLAGS) -c -o $@ $<
lens_light: $(LL_OBJ)
	$(GPP) $(LL_FLAGS) -I $(INC_DIR) -o $(BIN_DIR)/lens_light $(LL_OBJ) $(LL_LIBS)
lens_light_clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*
############## END: vkl_lens_light
