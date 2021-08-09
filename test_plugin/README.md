Test the Plugin (Command Line Mode)
===================================

Some sample files have been included in the ```test_plugin``` directory of this
repository.

To test the loading of a sample MCNP geometry and materials:

- Download or clone the repository
- Tthen navigate to the directory with ```cd Cubit-plugin/test_plugin```
- Then run the ```coreform_cubit``` command
- In the Cubit command line, type ```import mcnp test.i```
- Now save the geometry in ACIS format (.sat file) by entering ```export acis test.sat overwrite attributes_on``` in the Cubit command line.
- Finally enter ```export dagmc "test.h5m" faceting_tolerance 1e-3 make_watertight verbose``` into the cubit command line which will facet the geometry and save it in a format that can be used by DAGMC (.h5m file).

The same test can also be carried out from the regular command line without
loading the Cubit GUI by running journals (.jou files).

- In the regular commmand line enter ```coreform_cubit -batch -nographics -nojournal test.jou```
