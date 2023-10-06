General Capabilities (general)
    A) - Generate microlensing light curves on the fly for 2 bands
  * B) - Generate SN microlensing light curves on the fly for 2 bands
  * C) - Output cutouts
    D) - Just a static mock: SIE+shear, analytic source

No time variability (static) - main comparison general/D
    A) - Use a SPEMD instead of an SIE
    B) - Replace gaussian with equivalent Sersic sources
    C) - Input given in the COOLEST format
    D) - Use a pixelated source
    E) - Use a SIE and add potential perturbations
    F) - Use 2 lenses – one main and one satellite

AGN source (quasar) - main comparison general/A
    A) - Run with given AGN microlensing light curves
    B) - Change from mass-to-light ratio to a specific compact mass profile.
    C) - Add a non-microlensed variability component for 1 band
    D) - Use a custom accretion disc profile (still using the ‘moving fixed source’ model)
 $* E) - ‘variable moving source’ using an accretion disc ‘movie’ as input

Supernova source (supernova) - main comparison general/B
  * A) - Run with given SN microlensing light curves 
 $* B) - ‘variable moving source’ using a Supernova movie as input


* Not implemented yet
$ Computationally demanding
