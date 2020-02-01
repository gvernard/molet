GPP = g++


############## START: Fproject
CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  =  -ljsoncpp -lvkl -lgfortran -lCCfits -lcfitsio -lgmp -lCGAL
EXT_LIB_DIR = common/vkl_lib/lib

ROOT_DIR = lensed_point_source/vkl_point_source
EXT_INC_DIR = common/vkl_lib/inc
SRC_DIR = $(ROOT_DIR)/src
INC_DIR = $(ROOT_DIR)/inc
BIN_DIR = $(ROOT_DIR)/bin
OBJ_DIR = $(ROOT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))


DEPS = polygons.hpp contour-tracing.hpp simplify_caustics.hpp pointImage.hpp
OBJ  = polygons.o   contour-tracing.o   simplify_caustics.o   point_source.o
FULL_DEPS = $(patsubst %,$(INC_DIR)/%,$(DEPS)) #Pad names with dir
FULL_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir
#$(info $$OBJ is [${FULL_DEPS}])

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(FULL_DEPS)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -I $(EXT_INC_DIR) -c -o $@ $<

point_source: $(FULL_OBJ)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -I $(EXT_INC_DIR) -o $(BIN_DIR)/point_source $(FULL_OBJ) $(CPP_LIBS) -L $(EXT_LIB_DIR) -Wl,-rpath,$(EXT_LIB_DIR)
point_source_clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*
############## END: Fproject
