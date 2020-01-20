######## lenses:
A list holding all the different lenses.
There could be more than one lenses at the same or different redshift.
Each lens properties are divided into three groups: main_mass_model, substructure, and light_profile.

######## main_mass_model:
The type of the main (or smooth) mass model for the lens.
It could be analytical (SIE, SPEMD, etc), or from numerical simulations (e.g. EAGLE).
Each of these choices will result in an additional set of its own unique parameters.

######## substructure:
Superimposed substructure for each lens (at its defined redshift).
While the effects of individual perturbing analytical haloes can be considered as additional main lenses, this is not the case for extended perturbations.
These could include Gaussian Random Fields, or higher order moments in the lens mass distribution, such as a disc, spiral arms, etc.

######## light_profile:
The spatial distribution of source brightness, analytical or numerical.
One distribution per wavelength range.



######## source:
There can be only one source, which is basically just a 'light_profile'.
In case of multiple lens planes, the light profile of the first lens is sufficient to act as a source (together with the lensed features of the source due to this same lens).

######## point_source:
Specific quasar or SN properties for the source.




######## microlensing:
Single-plane microlensing for the moment.
