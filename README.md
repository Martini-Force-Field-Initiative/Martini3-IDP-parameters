# Supporting scripts and data for "Martini3-IDP: Improved Martini 3 Force Field for Disordered Proteins"

## DEPRECATION UPDATE JULY 2025

The contents of this repository has been fully integrated into both [Vermouth-Martinize](https://github.com/marrink-lab/vermouth-martinize)
and [Polyply](https://github.com/marrink-lab/polyply_1.0/). The tutorials detailed in [force_field/readme.md](force_field/readme.md) are deprecated.
In both programs, the force field described in the files here and the associated publication can be specified using `-lib martini3idp` (polyply) or `-ff martini3idp` (martinize2).
If you encounter problems with either implementation, open an issue in the respective repository.

In the case of martinize2, the program now handles disordered regions using the `-id-regions` and `-idr-tune` flags. No customised version
of Vermouth-Martinize is necessary to annotate these regions as previously described, and is _strongly_ discouraged.

### Reference
Wang, L., Brasnett, C., Borges-Araújo, L. et al. Martini3-IDP: improved Martini 3 force field for disordered proteins. *Nat Commun* 16, 2874 (2025). 
https://doi.org/10.1038/s41467-025-58199-2


### Martini3-IDP force field parameters and all related data/scripts.


- structure data for most simulations from the paper are available [here](https://zenodo.org/records/14608855)

- Martinize2 and Polyply format force field files are available in the [force_field](force_field) folder

- User tutorials are available [here](force_field/readme.md)

- Detailed processed data and python scripts of each testing are provided under corresponding folder
