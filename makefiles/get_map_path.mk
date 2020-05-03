.DEFAULT_GOAL := get_map_path

GPP = g++
CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  = -lgerlumph

ROOT_DIR = variability/extrinsic/get_map_path
SRC_DIR = $(ROOT_DIR)/src
BIN_DIR = $(ROOT_DIR)/bin
OBJ_DIR = $(ROOT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))

OBJ = get_map_path.o
FULL_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(GPP) $(CPP_FLAGS) -c -o $@ $<

get_map_path: $(FULL_OBJ)
	$(GPP) $(CPP_FLAGS) -o $(BIN_DIR)/get_map_path $(FULL_OBJ) $(CPP_LIBS)
clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*
