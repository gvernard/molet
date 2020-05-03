# LIBRARY: VKL_INSTRUMENTS
#======================================================
mkfile_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
vkl_instruments:
	make -f makefiles/vkl_instrument_modules.mk INSTRUMENT_PATH=$(mkfile_path)instrument_modules/ instrument_modules
vkl_instruments_clean:
	make -f makefiles/vkl_instrument_modules.mk clean


# ANGULAR_DIAMETER_DISTANCES
#======================================================
angular_diameter_distances:
	make -f makefiles/angular_diameter_distances.mk angular_diameter_distances
angular_diameter_distances_clean:
	make -f makefiles/angular_diameter_distances.mk clean

# VKL_FPROJECT
#======================================================
vkl_fproject:
	make -f makefiles/vkl_fproject.mk fproject
vkl_fproject_clean:
	make -f makefiles/vkl_fproject.mk clean

# VKL_QUASAR
#======================================================
vkl_point_source:
	make -f makefiles/vkl_point_source.mk point_source
vkl_point_source_clean:
	make -f makefiles/vkl_point_source.mk clean

# VKL_LLM
#======================================================
vkl_llm:
	make -f makefiles/vkl_llm.mk lens_light_mass
vkl_llm_clean:
	make -f makefiles/vkl_llm.mk clean

# GERLUMPH_MOVING_SOURCE
#======================================================
gerlumph_moving_source:
	make -f makefiles/gerlumph_moving_source.mk gerlumph_moving_source
gerlumph_moving_source_clean:
	make -f makefiles/gerlumph_moving_source.mk clean

# MATCH TO GERLUMPH MAPS
#======================================================
match_to_gerlumph:
	make -f makefiles/match_to_gerlumph.mk match_to_gerlumph
match_to_gerlumph_clean:
	make -f makefiles/match_to_gerlumph.mk clean

# GET_MAP_PATH
#======================================================
get_map_path:
	make -f makefiles/get_map_path.mk get_map_path
get_map_path_clean:
	make -f makefiles/get_map_path.mk clean

# INSTRUMENT_COMBINED_LIGHT
#======================================================
combined:
	make -f makefiles/combined_light.mk combined_light
combined_clean:
	make -f makefiles/combined_light.mk clean




ALL_DEPS := vkl_instruments
ALL_DEPS += angular_diameter_distances
ALL_DEPS += vkl_fproject
ALL_DEPS += vkl_point_source
ALL_DEPS += vkl_llm
ALL_DEPS += get_map_path
ALL_DEPS += match_to_gerlumph
ALL_DEPS += gerlumph_moving_source
ALL_DEPS += combined

CLEAN_DEPS = $(patsubst %,%_clean,$(ALL_DEPS))
#$(info $$OBJ is [${CLEAN_DEPS}])

all: $(ALL_DEPS)
clean: $(CLEAN_DEPS)
