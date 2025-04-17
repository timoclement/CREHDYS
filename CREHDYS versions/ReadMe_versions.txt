In this folder you can find the different versions of the Fortran code for the CREHDYS model.

------------------------------ CREHDYS V1 ------------------------------------------------------------------------------------------------------------------------------------------------------------------

- CREHDYS v1.0 is the original CREHDYS code as left by Eric Laloy on UCLouvain backups (Alfresco and physical hard drive) 
	- Running/compiling this code on Microsoft Visual Studio (MVS) with the Intel fortran compiler (ifort), using the input example files provided by Laloy, does not work properly and does not reproduce Laloy's original outputs.
There are convergence issues for the kinematic wave and NaN values in erosion outputs, among other problems.
	- In addition, there are several miscalculations in the code :
	 	- In the water content (theta) adjustement (daily and within events), the unit of control volume df is mistaken: it is treaed as [mm] while it is actually in [cm]
		- In "percol" and "drain" functions, there are unit errors: dthr should be 24h in percol (not dthr/3600), and dr is considered in [mm/h] while it's actually in [mm/dthr]. An unexpected "*thets" appears in the travel time calculation
		- A "*10" factor mistakenly inflates the settling velocity calculation (L354), making it 10 times too high
		- In the unit stream power equation, S is meant to be the SINE OF the slope, not just the slope itself
		- In the "etpr2" subroutine, "max1" is used but returns an integer instead of a real. Similarly, "min1" is incorrectly used elsewhere.
		- For long events, daily percolation is calculated again even though drainage was already accounted for during the event
		- There are error in evapotranspiration unit calculations (subroutine "etp")
		- If multiple cells have the maximum flow accumulation value (e.g. last row of a rectangular plot), output variables (e.g., total runoff and erosion) are computed only for the "last" outlet (bottom right cell), underestimating total fluxes.
		- Manning's n for channel flow (wheel tracks) is ignored - only overland flow Manning values are applied to all cells, including wheel tracks
		- The "etpr2" subroutine is executed before "edisp(m)(fc)" calculation, causing time/space shifts in evapotranspiration
		- In the ponding depth "s(m)" calculation, "filt(m)" is used instead of filtreal, which reflects actual infiltration rate as function of net rainfall and preceding ponding depth
		- For splash erosion on covered soil, full kinetic energy "ekbelcov2" is used, but it should be corrected for leaf interception (i.e., use ekbelcovsplash)

- CREHDYS v1.1 is the v1.0 original CREHDYS model, with some adjustements to allow successful compilation and execution in MVS using ifort. It produces outputs matching Laloy's examples.
	- qmin in the call to the "iter" (errq2) subroutine (used for Newton-Raphson iterations to solve kinematic wave outflow) is set to 0 instead of a small negative value, fixing convergence issues.
	- A condition checking for non-zero flow values was added to the sediment output calculations (in the explicit finite difference scheme), avoiding division by zero.
	- Recommended compiler options (Debug and Release):
		- Fortran>Optimization>Inline Function Expansion = Disable
		- Fortran>Language>Fixed Form Line Length = 132 Columns
		- Fortran>Data>Default KINDs: Integer/Real/Double Precision = 4/4/8, respectively
		- Fortran>Local Variable Storage = All Variable SAVE
		- Fortran>Floating-Point Exception Handling = "Underflow gives 0.0; Other exceptions produce NaN"
		- Fortran>Floating-Point Speculation = Fast
		- Fortran>Run-time>Check Array and String Bounds = No
		- Fortran>Run-time>Check Stack Frame = No
	- Additional Debug-only options:
		- Fortran>General>Debug Information Format = Full ('None' in "Release")
		- Fortran>General>Optimization = Disable ('Maximize speed' in "Release")
		- Fortran>Diagnostics>Check Routine Interfaces = Yes
		- Fortran>Run-time>Generate Traceback Information = Yes
		- Fortran>Run-time>Check Stack Frame = Yes


------------------------------ CREHDYS V2 ------------------------------------------------------------------------------------------------------------------------------------------------------------------

- CREHDYS v2.0 is the CREHDYS v1.1 model, with corrections of all major errors mentioned above and other minor adaptations.
	- Major corrections :
		- df units converted from [cm] to [mm] for proper water content adjustement
		- variables and outputs of "drain" and "percol" function/subroutine corrected (units, *thets, dthr, etc.)
		- Corrected sediment settling velocity formula (L354)
		- Used sine of slope in unit stream power calculation (L1727)
		- Replaced min1/max1 with standard min/max in "etpr2" subroutine and elsewhere
		- Daily percolation is skipped if a rain event already triggered drainage that day. However, rain events in weather.inp(ut) must now last 24h, with 0 mm/h intensity during dry periods to allow drainage calculation.
		- Evapotranspiration unit errors fixed
		- "totrun" and "toteros" now aggregate all outlet contributions, not just the last cell with the maximum flow acc
		- Manning's n values now vary depending on cell type (overland flow OR wheel track)
		- "etpr2" subroutine is now called after edisp(m)(fc) calculation

	- Minor adptations :
		- New outputs files are generated to follow daily dynamic of crop and soil systems :
			- "crop.out": (inter)crop/maize growth and cover
			- "soil.out" for soil hydrological parameters and water fluxes (infiltration, drainage, evapotranspiration, water content)
		- The "eps" tolerance in "gamptun3" (Newton-Raphson iterations for cumulative infiltration) was increased from 1.0e-6 to 2.0e-6 to avoid infinite while loops caused by limited precision


- CREHDYS v2.1 is the CREHDYS v2.0 model, with adaptations for the "Intell'eau" project.
	- New features for Intell'eau :
		- Crop growth model ajusted to simulate simultaneous maize growth and cover crop decay (first season only). Intercrop residues, vegetation cover and LAI calculations are updated: dead intercrop residues are ignored in LAI (interception and transpiration) but still protect the soil from sealing
		- Percolation is now calculated before AND after each rain event:
			- This speeds up simulations (no need to compute full-day zero rainfall for drainage only)
			- This should provide more realistic infiltration rate, based on initial soil moisture.
			- NB : Drainage is still also calculated during the rainfall event
	- Other minor corrections :
		- "filtreal" replaces "filt(m)" for ponding depth calculation
		- Weighting factor in the implicit(-explicit) scheme "Omega" was increased from 0.75 to 0.80 to reduce numerical oscillations of outflow "q22" due to kinematic shocks occuring at abrupt changes in rainfall intensity. This choice is based on recommendation from the KINEROS documentation.
		- Totl sand (coarse sand "sa(i)" + very fine sand "vfs") in the pedotransfer functions (max. evaporative depth and soil transmissivity) 
		- In "etpr2" subroutine, evaporation follows original Ritchie (1972)'s rule: back to stage 1 when cumulative rainfall under stage 2 > cumulative evaporated water under stage 2. (previously based on water content > field capacity, leading to long stage 2 with low evaporation values even on moist soil.
		- Condensed code and fiwed issues in percolation, evaporation and transpiration routines.
		- Corrected KE used for splash erosion on covered soil (rain falling of maize leaves "ekbelcovsplash")

- CREHDYS v2.2 is the CREHDYS v2.1 model, but modified by Thomas de Maet for successful compilation with gfortran on Linux (e.g., use on CISM supercomputer clusters).

------------------------------ CREHDYS V3 ------------------------------------------------------------------------------------------------------------------------------------------------------------------

- CREHDYS v3.0 is the CREHDYS v2.2 model, with a correction to modeling of the splash erosion process :
	- Splash detachment is now calculated separately for dry (water depth=0mm) and ponded areas, as in the openLISEM model.

------------------------------ CREHDYS V4 ------------------------------------------------------------------------------------------------------------------------------------------------------------------

- CREHDYS v4.0 is the CREHDYS v3.0 model, adapted for micro-basin tillage modeling :
	- "MB_storage" variable added to represent the surface storage capacity in the micro-basins, in addition to the surface storage in the primary microtopography which is calculated from random roughness
	- Connectivity threshold "connectresh" added: runoff is triggered when water height exceeds connectresh*total surface storage
	- The equation for surface storage in the microtopography was updated to that from Kamphorst et al. (2000), which is used in openLISEM
	- Equation for surface water storage (L1477) corrected to include runoff from upstream even when there is no runoff in the cell, via "+ hindt"
