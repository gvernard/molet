.DEFAULT_GOAL := angular_diameter_distances

GPP = g++


CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  = -ljsoncpp

ROOT_DIR = cosmology/angular_diameter_distances
INC_DIR = $(ROOT_DIR)/inc
SRC_DIR = $(ROOT_DIR)/src
BIN_DIR = $(ROOT_DIR)/bin
OBJ_DIR = $(ROOT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))


DEPS = auxiliary_functions.hpp
OBJ  = auxiliary_functions.o   angular_diameter_distances.o
FULL_DEPS = $(patsubst %,$(INC_DIR)/%,$(DEPS)) #Pad names with dir
FULL_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir
#$(info $$OBJ is [${FULL_DEPS}])

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(FULL_DEPS)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -c -o $@ $<

angular_diameter_distances: $(FULL_OBJ)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -o $(BIN_DIR)/angular_diameter_distances $(FULL_OBJ) $(CPP_LIBS)
clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*

