# VKL_LIB
vkl_lib:
	make -f makefiles/vkl_lib.mk vkl_lib
vkl_lib_clean:
	make -f makefiles/vkl_lib.mk vkl_lib_clean

# VKL_FPROJECT
vkl_fproject: vkl_lib
	make -f makefiles/vkl_fproject.mk fproject
vkl_fproject_clean:
	make -f makefiles/vkl_fproject.mk fproject_clean

# VKL_QUASAR
vkl_point_source: vkl_lib
	make -f makefiles/vkl_point_source.mk point_source
vkl_point_source_clean:
	make -f makefiles/vkl_point_source.mk point_source_clean

# VKL_LLM
vkl_llm: vkl_lib
	make -f makefiles/vkl_llm.mk lens_light_mass
vkl_llm_clean:
	make -f makefiles/vkl_llm.mk lens_light_mass_clean


# GERLUMPH_LIB
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))
gerlumph_lib:
	make -f makefiles/gerlumph_lib.mk cpu MAP_PATH=$(current_dir)/data/gerlumph_mirror/
gerlumph_lib_clean:
	make -f makefiles/gerlumph_lib.mk clean

# GERLUMPH_MOVING_SOURCE
gerlumph_moving_source: gerlumph_lib
	make -f makefiles/gerlumph_moving_source.mk gerlumph_moving_source
gerlumph_moving_source_clean:
	make -f makefiles/gerlumph_moving_source.mk gerlumph_moving_source_clean




all: vkl_fproject vkl_point_source vkl_llm
clean: vkl_lib_clean vkl_fproject_clean vkl_point_source_clean vkl_llm_clean
