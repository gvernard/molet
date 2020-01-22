# VKL_FPROJECT
vkl_fproject:
	make -f makefiles/vkl_fproject.mk fproject
vkl_fproject_clean:
	make -f makefiles/vkl_fproject.mk fproject_clean

# VKL_LENS_LIGHT
vkl_lens_light:
	make -f makefiles/vkl_lens_light.mk lens_light
vkl_lens_light_clean:
	make -f makefiles/vkl_lens_light.mk lens_light_clean




all: vkl_fproject vkl_lens_light
clean: vkl_fproject_clean vkl_lens_light_clean
