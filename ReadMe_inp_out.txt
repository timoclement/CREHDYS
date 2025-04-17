Details of .inp and .out files

1. gen.inp
Plot name
Simulation duration [# of days]
Computation time step (within events)
Number of soils
Total porosity [m³/m³soil], field capacity and wilting point (as fraction of total porosity), Ksat of subsoil [mm/h], depth of topsoil layer (control volume) [cm]
Clay, sand, silt, coarse fragments, very fine sand [%]
Random roughness in wheel track [mm], Manning's n in wheel track [s.m^(-1/3)]
Soil albedo [/]
Cell size [cm]
Number of intercrop and maize cropping season periods
Starting day(s) of intercrop period(s)
Intercrop destruction and burial dates
Starting day(s) of maize period(s)
Number of wheel tracks in intercrop and maize periods
Initial random roughness [mm] for intercrop period(s)
Initial random roughness for maize period(s)
Initial Ksat [mm/h] for intercropping and maize periods
Dummy indicating if RR and Ksat are varying (=1) or not (=0) during each intercropping period
Dummy indicating if RR and Ksat are varying (=1) or not (=0) during each maize period
Cover extinction coefficient [m-2/m²] in eq LAI=f(cov) for intercrop and maize periods
Growth rate coefficients for intercrop and maize
Coefficient for cover crop cover=f(biom)

2. plot.inp
line1: number of grid elements, number of rows, number of columns
line2 -> end: element index, row, column, slope (%), flow direction, flow accumulation, soil index, crop index, channel index (0 = overland flow cell, 1 = wheel track).
-9999 means element is outside of the basin limits

3. weather. inp
line1: day of the simulation, daily surface air temperature (°C), daily surface soil temperature (not used), daily solar radiation (MJ/m2), dummy (1 = a rainfall event will occur this day, 0 = no rainfall event), time of event start [hour since the day started].
line2: if rainfall event occurs this day: day of the simulation, time in the event (min), rainfall rate (mm/hr), dummy (not used), dummy (1=end of the event, 0= event not ended yet), dummy (not used).
If no rainfall, same as line1.

4. theta.inp
theta: Initial soil moisture of grid elements (m³ water/m³ soil) with respect to their index,
-9999 means element outside of the basin limits.

5. paramhyd.inp
line1:
kso2: final overland flow cell saturated hydraulic conductivity (mm/hr),
Kcho2: saturated hydraulic conductivity of wheel tracks during the intercropping period (mm/hr),
kchm2: saturated hydraulic conductivity of wheel tracks during the maize period (mm/hr),
hfront: Green-Ampt soil matric potential at the wetting front (mm),
c_kso: Soil stability factor for the time-evolution of the overland flow cell saturated hydraulic conductivity (m2/J),
c_rr: Soil stability factor for the time-evolution of the overland flow cell random roughness (m2/J),
mann: Manning's n coefficient for overland flow cells (s. m^-1/3)
MB_stor: Surface storage in microdepressions due to micro-basin tillage [mm]
prop_nodep: Proportion of surface that is not under microdepression [/]
connectresh: Connectivity threshold = proportion of maximum total surface storage (primary microtopography + micro-basins) at which runoff is initiated [/]

6. parameros.inp
line1:
As: Number of drops of the aggregate stability test of Low (1954) (-),
d50: Median soil particle size (µm)
Coh: Soil cohesion of overland flow cells measured by a torvane (kPa)
Cohw: Soil cohesion of wheel track cells measured by a torvane (kPa)

7. daily.out
line1-> ...
day, daily rainfall (mm), daily runoff (total for all outlets) (mm), daily erosion (total for all outlets) (t/ha), theta at outlet at end of day
3 last lines:
total rainfall, runoff and soil losses over the simulation

8. dailyRE.out
day, daily rainfall (mm), daily runoff (total for all outlets) (mm), daily erosion (total for all outlets) (t/ha)

9. eventR.out
day, rainfall rate (mm/hr), net rainfall rate (after interception) [mm/h], computation time (min),
(potential) infiltration rate [mm/h], mean (for all outlets) ponding depth [mm], mean surface retention (connectresh*DSmax) [mm],
total outflow from all outlets relative to whole plot area (mm/hr), total outflow from all outlets (m3/s), mean runoff height of all outlets (mm), 
theta at last outlet

10. eventE.out
day, rainfall rate (mm/hr), time (min), total sediment flux from all outlets (kg/s), 
Transport capacity at the last outlet (g/l), mean sediment concentration in the runoff water from all outlets (g/l),
deposition at the last outlet (kg/s), flow detachment at the last outlet (kg/s), splash detachment at the last outlet (kg/s),
net erosion at the last outlet [kg/s].

11. crop.out
day, current cropping period (maize or intercrop), daily average temperature [°C],
intercrop (residue) biomass [g/m²], cover of intercrop (residues) [%], cover of maize [%],
LAI of intercrop, LAI of maize

12. soil.out
day, daily rainfall [mm], (dynamic) Ksat at last outlet [mm/h], RR at last outlet [mm],
effectively infiltrated water at last outlet [mm/day],
effective drainage rate at last outlet [mm/day], 
actual evaporation rate at last outlet [mm/day], actual transpiration rate at last outlet [mm/day],
total effective water loss from evapotranspiration and percolation [mm/day],
soil water content [/]