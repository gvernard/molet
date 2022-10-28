Here is the summary of all the examples included:

"paper_data" -  All the cases presented in the MOLET paper (Vernardos 2021).
"standard"   -  An example using a .coolest file as input.
"test_0"     -  A basic test for an AGN point source, with given intrinsic and microlensing variability light curves.
"test_1"     -  Same as 'test_0', but it includes an additional unmicrolensed intrinsic signal from the source (e.g. coming from a large area like the BELR).
"test_2"     -  Same as 'test_0', but now the microlensing light curves will be generated from magnification maps on the fly, using the 'moving_fixed_source' model. If the maps are missing, a download link is provided automatically.
"test_3"     -  A selected microlensing light curve trajectory from 'test_0' (index: 20), for which cutouts are also calculated.
"test_4"     -  Same as 'test_0' but now the point source is replaced by an 'expanding_source' and the cadence of the observations is changed.
"test_5"     -  Same as 'test_0' but now the point source is removed and there are gridded lens potential perturbations present.
"test_6"     -  Same as 'test 2' but instead of the 'moving_fixed_source' variability model we have the 'moving_variable_source' model.
"test_7"     -  Same as 'test 0' but instead of the an 'SIE' we use a 'SPEMD' model for the lens mass.
