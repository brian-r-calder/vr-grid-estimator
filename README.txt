Center for Coastal and Ocean Mapping & NOAA/UNH Joint Hydrographic Center
University of New Hampshire
Durham, NH 03820

May 23, 2017

This example code distribution accompanies the paper "Computationally Efficient Variable Resolution Depth Estimation", B.R. Calder and G. Rice, published in Computers and Geosciences.  The code demonstrates the ideas expressed in the paper, implementing a 1D version of the algorithm in order to keep the level of complexity to a minimum, which hopefully makes it easier to understand in addition.  The code is implemented in MATLAB, but could be converted to other languages relatively simply.  Note that because of the implementation method, the code is not expected to be run-time efficient, although efforts have been made to keep the runtime within reasonable bounds on modern desktop machines.

The code is organised in a set of class objects that represent the core ideas within the paper, a series of helper scripts (e.g., to display the results of some stages of the computation), and an example estimation script which ties together all of the components to demonstrate a simple example of what an estimation sequence would look like in practice.  Summary details of the files are given below.

This code is copyright of the University of New Hampshire as indicated in the source files, and is released under version two of the GNU General Public License.  Please see LICENSE.txt for the full text of the license.

File Description
------------------

Example Estimation Script
--------------------------------

In order to demonstrate the usual method for conducting estimates with these algorithms, the example_estimation.m script ties together all of the components of the code, from generating a simulated bathymetric profile to displaying the final estimates in comparison with the bathymetric model.  It also demonstrates how to call the other objects provided in the example code.

The "Simulator Parameters" cell provides the parameters for the simulated data generation, and is the most likely area for modifications during initial testing.  The simulation allows depth to be positive up or down, but to allow for simpler plotting of results, the script uses positive up depths (hence the (numerically) minimum depth is more negative than the (numerically) maximum depth).  Other plausible modifications are the roughness of the bathymetry (which must be in the range between 1.0 and 2.0), the ship speed, and sounder beamwidth.  Note that the beamwidth is set in radians; depending on your version of MATLAB and the toolboxes installed, you may have to provide the deg2rad() functionality.

The "Estimator Parameters" cell provides the configuration for the TileManager and estimators.  Most of the parameters are either self-explanatory, or discussed in the accompanying paper.  Note that the code uses a simple numerical value to select which underlying estimator is being used (MeanEstimator = 1, WeightedMeanEstimator = 2).

The remaining sections of the script should not need much, if any, modification, since they use only core MATLAB functionality and the classes provided in the example code.  In order to allow for detailed inspection of the algorithm's operation, the script supports three external variables which adjust some of the behaviours.  First, setting 'record' true in the global namespace before running will cause the script to keep a record of the observations added to the estimator.  Second, setting 'preserve' true will cause the script to keep a copy of all of the EstimatorTile objects used during the first pass of the estimator (which are normally replaced with RefinedTile objects prior to the second pass starting).  Finally, setting 'rerun' true after a completed run of the estimator will cause the script to use the preserved PRNG state from the first run, rather than re-seeding the PRNG.  This causes the code to repeat the same sequence of observations, therefore allowing the same computation to be run again, rather than generating a new bathymetry sample.

When complete, the code will generate a total of six figures: three for the first-pass observation density and sample spacing estimates, one for the high-resolution refined estimates, and two comparing the estimated depth with the original model.

Class Objects
----------------

Estimator Hierarchy: Estimator.m, MeanEstimator.m, and WeightedMeanEstimator.m

Estimator.m provides the Estimator abstract base class for anything that provides estimation services at a single point in the computational domain; MeanEstimator and WeightedMeanEstimator are derived from this and represent simple sample mean, and uncertainty-weighted sample mean estimation, respectively.  All estimators also compute a sample estimate of the standard deviation of the observations presented against the eventual mean value.  In order to avoid having to copy the estimator object as it is updated with new observations, Estimator is a handle object in MATLAB-speak (effectively, passed by reference).

The Estimator class is provided so that interested readers could substitute their own estimator without modifying the remainder of the code.  Anything that supports the abstract interface in Estimator should be able to be substituted directly; edit the Estimator.Create() method in Estimator.m to recognise the new estimator's ID number, and instantiate it appropriately.

Tile Hierarchy: Tile.m, EstimatorTile.m, RefinedTile.m, and DummyTile.m

Tile.m provides the Tile abstract base class to represent the central concept of a contiguous group of low-resolution estimation cells somewhere in the computational domain, which are handled as a group.  The size of the cells, and the cell-count per tile, are user parameters.  Splitting the computational domain into tiles allows TileManager (below) to turn on resources only where data is observed.

Sub-classes of Tile provide for different contents in the estimation cells, and therefore for first-pass and second-pass processing.  The EstimatorTile sub-class provides first-pass estimation of observation density, and therefore predicted high-resolution refined grid sample spacing for each estimation cell in the tile.  The RefinedTile sub-class provides second-pass estimation of depth (in this example), based on the estimates of refined grid sample spacing provided by the corresponding EstimatorTile.  The design here is that a RefinedTile is constructed from, and then effectively replaces, an EstimatorTile in the same physical location.

DummyTile is used to provide a null-constructor compatible Tile sub-class; this allows for Tile sub-classes to be held in heterogeneous arrays (MATLAB-speak for a polymorphic array of references); Tile has the appropriate decorations to allow for this, and generates DummyTile objects when asked for a default object.

Georeferencing for a Tile is handled locally (i.e., the tile has no sense of its absolute location in the computational domain, assuming that the code passing in the observations makes the appropriate offsets).  Georeferencing for an individual estimation cell within the tile is at the centre of the estimation cell (i.e., the left-most estimation cell is a half cell-width from the zero point).  In RefinedTile estimation cells, the refined grid is always centred in the estimation cell.

For efficiency, Tile is a MATLAB handle class, and is therefore effectively passed by reference, and can be updated in place.

Overall Computation: TileManager.m

TileManager.m contains the TileManager class, which provides estimation services over a unified computational domain; it is the primary object that is instantiated to manage the various Tile sub-classes required to do a computation using these methods.  TileManager assumes a uniform coordinate frame in the input observations, where the ordinates must be positive and increasing from left to right across the computational domain.  Otherwise, where the observations occur is irrelevant, and the size of the domain over which the observations occur is limited only by memory available.  (Note that this version of the algorithms, for the sake of simplicity and explicability, do not implement the memory-mapped data structures used for the 2D production version of the code as described in the paper).

User-level code should mostly use TileManager as its interface to the algorithms, although the individual components (e.g., Tile sub-classes within the TileManager cache) can be accessed for further study.  The TileManager.Dimensions method provides a summary of the current state of the tile cache.

Note that TileManager maintains only one Tile sub-class object for each tile in the computational domain.   Hence, if a tile is refined, the EstimatorTile is replaced with a RefinedTile, and the initial estimates are deleted.  If observations are added to an area where no observations have been seen previously, an EstimatorTile is constructed for them by default.  Consequently, it is possible to have a mixture of EstimatorTile and RefinedTile objects in the cache at the same time, although in most cases of expected usage this should be rare.  Calls that reference first-pass (observation density) estimates only use EstimatorTile objects for their output.  Calls that reference generic "depth" estimates use the low-resolution auxiliary (fixed spacing) depth estimates constructed by EstimatorTile objects at the estimation cell resolution, and high-resolution (variable spacing) depth estimates constructed by RefinedTile objects.  In either case, the spacing of the estimates is also given in the output matrix, allowing the user-level code to determine which is which (if required).

Georeferencing for a tile is constructed by the TileManager (Tile sub-classes hold no internal sense of where they are in the computational domain: the TileManager must offset observation ordinates for them).  The georeferencing spot for a tile is at the left hand edge of the area covered.

Input Observations: Observation.m

Observation.m provides the Observation class, which encapsulates all of the required information about an observation, except its location.  Note that the uncertainty is recorded as a variance in order to avoid a number of square-root and square operations during computation.

Simulated Data: GenerateBathymetry.m and SimulateObservations.m

In order to provide for controlled test data that can be used to examine the behaviour of the algorithms, GenerateBathymetry.m provides the GenerateBathymetry object, which can sample random depth profiles using a power-spectrum fractal data generator.  User parameters control the degree of roughness, the dynamic range, and a filtering parameter to exclude very small wavelength features, if required.  The code arranges the generation of the fractal with over-sampling in order to avoid first-harmonic distortion, and aliasing.

SimulateObservations.m provides the SimulateObservations object in order to transform the bathymetry generated by GenerateBathymetry.Sample into a sequence of observations such as might be generated by a simple echosounder moving over the bathymetry at a constant speed-through-water in an area of constant geometric mean sound speed.  The simulation model is extremely simple, assuming that the echosounder generates only one ping in the water at a time, and sends the next ping as soon as the first is received.  The simulation assumes that the echosounder has fixed beamwidth, but the Observation object allows the beamwidth to be specified separately for each observation, e.g., to emulate different beams from a multibeam echosounder that are being used along-track.  More complex models are obviously possible.

Helper Scripts
-----------------

Tile-based Visualisation: summarise_esttile.m and summarise_reftile.m

These scripts generate a sequence of plots to display the current contents of an EstimatorTile and RefinedTile, respectively.  The code attempts to provide plausible decoration and formatting of the plots, but may need to be modified for specific applications.

TileManager Visualisation: display_density_estimates.m and display_refined_estimates.m

These scripts can be used to generate plots of the results extractable from a TileManager at the end of the first pass of an estimation, and the second pass, respectively.  The TileManager object returns the estimates extracted from the various tiles in cache (across a user-specified range, or the whole range if none is specified at input) as a unified matrix with estimate locations in the first column, and various estimates or statistics in the remaining columns.  Details of the contents of these matrices are provided in the documentation in the TileManager class file.
