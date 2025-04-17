# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 09:58:56 2024

@author: ticlement
"""

#%% CELL 1 : INPUT PREPARATION -------------------------------------------------------------------------------------------------------

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import csv
import pandas as pd
from datetime import datetime

current_wd = os.getcwd()
print("Current working directory :", current_wd)
os.chdir("C://Users//ticlement//Documents//MODELISATION//CREHDYS//CREHDYS Github") #Adapt path
current_wd = os.getcwd()
print("Current working directory :", current_wd)

# prep input GEN.INP #
plot_name="Ittre_2022"
simdur=124
dt=60 ; nsoil=1 #simdur [days], dt [s], nsoil (# of different soil types) can actually not be>1
tp=0.411 ; fc=0.609 ; wp=0.255 ; #total porosity [m³/m³ soil], field capacity and wilting point [as fraction of total porosity]
ksat_sub=20 ; df=20 #ksat for subsurface (drainage) [mm/h], df is the control volume, i.e., the depth of topsoil layer [cm]
cl=9 ; sa=3 ; st=36 ; wcf=0 ; vfs=52 # fractions of clay, sand, silt, coarse fragments and very fine sand [%]
rrch=5 ; mnch=0.03 ; #random roughness on channel (wheel tracks) [mm] and manning on channel [m^-1/3 / s]
albed=0.2 ; #soil albedo
dx=0.25 #cell size [m]
ninter=1 ; nmaize=1 #number of intercrop and maize periods. Not used.
startdayinter=[0] ; startdaymaize=[1] #dates of start of intercropping and maize periods [# of days since start of simulation]. Let zero if no period
dest=0 ; bur=0 #destruction timing (# of days since intercrop sowing) and burial timing (# of days since destruction) let 0 if no destruction/burial date
wheelinter=[0] ; wheelmaize=[0] #not used
rrini_inter=[9] ; rrini_maize=[9] #initial random roughness [mm] for intercrop and maize periods. Let 0 if no period
kinio=75 ; kinim=75 #initial Ksat for intercropping and maize period respectively [mm/h]
rrksvarint=[0] ; rrksvarmaize=[0] #dummy indicating if RR and Ksat are varying (=1) or constant (=0) during each intercropping and maize period
kextint=0.4 ; kextmaize=0.34 #cover extinction coefficient [m-2/m2] in eq. LAI=f(cov)
pcoverint=0.00008 ; kcovermaize=0.005 #growth rate coefficients for intercrop and maize
rbiomcover=0.0076 #coefficient for cover crop cover=f(biom)

from input_prep import prep_gen
prep_gen(current_wd,plot_name,simdur,dt,nsoil,tp,fc,wp,ksat_sub,df,cl,sa,st,wcf,vfs,rrch,mnch,albed,dx,ninter,nmaize,startdayinter,startdaymaize,dest,bur,wheelinter,wheelmaize,rrini_inter,rrini_maize,kinio,kinim,rrksvarint,rrksvarmaize,kextint,kextmaize,pcoverint,kcovermaize,rbiomcover)

# prep input PARAMHYD.INP #
kso2=1.91 ; kcho2=0.1 ; kchm2=0.1 #final Ksat on overland flow cells [mm/h]; Ksat [mm/h] of wheel tracks during intercropping and maize periods respectively
hfront=147 #Green-Ampt soil matric potential at wetting front [mm]
c_kso=0.002 ; c_rr=0.0001  #soil stability factor [m²/J] for temporal dymanics of Ksat and RR of overland flow cells
mann=0.016 #Manning'n coefficient for overland flow cells [s.m^-1/3]
MB_stor=0 #surface storage in microdepressions due to micro-basin tillage [mm]. Let 0 if conventionnal tillage
prop_nodep=1 #proportion of surface that is not under microdepression, i.e. where random roughness contribute to surface retention. Let 1 if conventionnal tillage
connectresh=1 #connectivity threshold: proportion of the total surface storage (primary microtopography + micro-basins) at which runoff is initiated. =0.5 according to Penuela-Fernandez (2015).

from input_prep import prep_paramhyd
prep_paramhyd(current_wd,kso2,kcho2,kchm2,hfront,c_kso,c_rr,mann,MB_stor,prop_nodep,connectresh)

# prep input PARAMEROS.INP #
as_=10 #number of drops of the aggregate stability test of Low (1954)
d50=60 #median soil particle diameter [µm]
coh=10 ; cohw=20 #soil cohesion [kPa] measured by a torvane, for overland flow cells and wheel tracks respectively

from input_prep import prep_parameros
prep_parameros(current_wd,as_,d50,coh,cohw)

# prep input PLOT.INP # #prepare a rectangular plot with homogeneous slope in y direction
xlength=20 ; ywidth=0.25 # width and length of a rectangular plot [m]. must be a multiple of dx
slope=13.9  #mean slope steepness [%]
yvert_wt=[] #y center position [m] of vertical wheel track(s) (e.g., 0.375). let empty list if no WT.
wdth_wt=float(0.25) #width of wheel tracks [m]. must be a multiple of dx

from input_prep import prep_rec_hom_plot
prep_rec_hom_plot(current_wd,dx,xlength,ywidth,slope,yvert_wt,wdth_wt)

# prep input THETA #
theta_ini=0.20 # initial soil moisture content [m³water/m³soil]

from input_prep import prep_theta
prep_theta(current_wd,dx,xlength,ywidth,theta_ini)

# prep input WEATHER #
from input_prep import prep_weather
prep_weather(current_wd,simdur,'Ittre 2022 daily weather.xlsx','Ittre 2022 event rain.xlsx')
 # Give as input of prep_weather the name of (xlsx) files containing daily weather (typical pameseb output with tsa, tsf (also tsa if tsf not available), plu, and rad (but a new rad column has to be created for conversion in MJ/m²)) and event rainfall data (typical rain gauge data with 2 columns : date, and rain).
 # See example for adequate order of columns
 
 
 
#%% CELL 2 : RUN CREHDYS MODEL -------------------------------------------------------------------------------------------------------
import subprocess

subprocess.run('CREHDYS_v4.0.exe', shell=True)  
dailyRE = pd.read_csv('dailyRE.out', sep='\s+', names=['day','rain','runoff','erosion'],header=None)



#%% CELL 3 : PLOT OUTPUTS -------------------------------------------------------------------------------------------------------
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import scipy.optimize

current_wd = os.getcwd()
print("Current working directory :", current_wd)
os.chdir("C://Users//ticlement//Documents//MODELISATION//CREHDYS//CREHDYS Github")
current_wd = os.getcwd()
print("Current working directory :", current_wd)

crop = pd.read_csv('crop.out', sep='\s+', names=['day','croperiod','temp','intbiom','covint','covmaize','LAIint','LAImaize'],header=None) # '\s+' as separator enables one single OR several spaces as possible separators
dailyRE = pd.read_csv('dailyRE.out', sep='\s+', names=['day','rain','runoff','erosion'],header=None)
soil = pd.read_csv('soil.out', sep='\s+', names=['day','rain','Ksat','RR','inf','drain','evap','transp','totwloss','theta'],header=None)
eventR = pd.read_csv('eventR.out', sep='\s+', names=['day','rainrate',',netrainrate','time','inf','s','dirvon','qmmhr','qm3s','hmm','theta'], header=None)

#add dates to dataframes
#starting day of simulation


# Plot CREHDYS outputs. /!\ these instructions have to be launched simultaneously (and not line by line) for proper plot display

#Plot soil.out
soil.plot(subplots=True, figsize=(10,20)) ; plt.legend(loc='best')

fig, axes = plt.subplots(3,1, figsize=(15,15))
plt.suptitle('Daily water fluxes in/out the soil', fontsize=20)
plt.rc
axes[0].bar(soil['day'],soil['rain'], width=0.5, label='rain [mm/d]')
axes[0].bar(soil['day'],-soil['inf'],color='darkturquoise', width=0.4, label='infiltration [mm/d]')
axes[0].bar(soil['day'],-soil['drain'],color='firebrick', width=0.4, label='drainage [mm/d]')
axes[0].scatter(dailyRE['day'],dailyRE['runoff'].replace(0,np.nan),marker='o', c='darkorange', label='runoff [mm/d]') #we do not want to see zero runoff values
axes[0].axhline(y=0, color='black', linestyle=(0,(5,5)), linewidth=0.8) #linestyle=(0,(5,5)) corresponds to a (loosely) dashed line
axes[0].set_ylabel('soil water fluxes [mm/day]', fontsize=15)
axes[0].legend(loc='lower left', fontsize=15)
axes[0].tick_params(labelsize=15)
axes0 = axes[0].twinx()
axes0.plot(soil['day'],soil['Ksat'],color='yellowgreen', label='Ksat [mm/h]')
axes0.plot(soil['day'],10*soil['RR'],color='violet', label='RR [mm*10]')
axes0.set_ylabel('Ksat [mm/h], RR [mm*10]', fontsize=15)
axes0.tick_params(labelsize=15)
axes0.legend(loc='upper right', fontsize=15)
fig.tight_layout()
fig.show()
axes[1].bar(soil['day'],soil['rain'], width=0.5, label='rain [mm/d]')
axes[1].set_ylabel('rain [mm/day]', fontsize=15)
axes[1].tick_params(labelsize=15)
axes1 = axes[1].twinx() #secondary y axis
axes1.plot(soil['day'],soil['theta'], color='saddlebrown', label='soil moisture [/]')
axes1.set_ylabel('soil moisture theta [/]', fontsize=15)
axes1.tick_params(labelsize=15)
axes[1].legend(loc='upper left', fontsize=15)
axes1.legend(loc='upper right', fontsize=15)
fig.tight_layout()
fig.show()
axes[2].plot(soil['day'],soil['evap'],color='chocolate', label='actual soil evaporation [mm/d]')
axes[2].plot(soil['day'],soil['transp'],color='forestgreen', label='actual plant transpiration [mm/d]')
axes[2].set_ylabel('potential evaporation & transpiration [mm/day]', fontsize=15)
axes[2].tick_params(labelsize=15)
axes[2].legend(loc='best', fontsize=15)
axes[2].set_xlabel('Day since start of simulation', fontsize=15) # Set x label only on last subplot because subplots share same x axis
fig.tight_layout()
fig.show()

# Plot crop.out
fig, axes = plt.subplots(2,1, figsize=(15,15))
plt.suptitle('(Inter)crop growth and cover', fontsize=20)
plt.rc
axes[0].scatter(crop['day'],crop['temp'], color='darkred', label='daily avg temperature [°C]')
axes[0].plot(crop['day'],crop['LAImaize'],color='gold', linewidth=1.6, label='maize LAI [/]')
axes[0].set_ylabel('daily temperature [°C] and maize LAI [/]', fontsize=15)
axes[0].legend(loc='upper left', fontsize=15)
axes[0].tick_params(labelsize=15)
axes0=axes[0].twinx() # secondary y axis
axes0.plot(crop['day'],crop['covmaize'],color='forestgreen', linewidth=1.6, label='maize cover [%]')
axes0.set_ylabel('maize cover [%]', fontsize=15)
axes0.tick_params(labelsize=15)
axes0.legend(loc='upper right', fontsize=15)
fig.tight_layout()
fig.show()
axes[1].scatter(crop['day'],crop['temp'], color='darkred', label='daily avg temperature [°C]')
axes[1].plot(crop['day'],crop['intbiom']/100,color='blueviolet', linewidth=1.6, label='intercrop biomass [t/ha]')
axes[1].plot(crop['day'],crop['LAIint'],color='gold', linewidth=1.6, label='intercrop LAI [/]')
axes[1].tick_params(labelsize=15)
axes[1].set_ylabel('daily temperature [°C], intercrop biomass [t/ha] and intercrop LAI [/]', fontsize=15)
axes[1].legend(loc='upper left', fontsize=15)
axes[1].set_xlabel('Day since start of simulation', fontsize=15) # Set x label only on last subplot because subplots share same x axis
axes1=axes[1].twinx() # secondary y axis
axes1.plot(crop['day'],crop['covint'],color='forestgreen', linewidth=1.6, label='intercrop cover [%]')
axes1.set_ylabel('intercrop cover [%]', fontsize=15)
axes1.tick_params(labelsize=15)
axes1.legend(loc='upper right', fontsize=15)
fig.tight_layout()
fig.show()

#Plot runoff and erosion outputs
fig, ax1 = plt.subplots(figsize=(14,10))
plt.suptitle('Daily surface fluxes', fontsize=20)
ax1.bar(dailyRE['day'],dailyRE['rain'], width=0.5, label='rain [mm/d]')
ax1.set_ylabel('daily rain [mm/day]', fontsize=15)
ax1.legend(loc='upper left', fontsize=15)
ax1.tick_params(labelsize=15)
ax2=ax1.twinx() # secondary y axis
ax2.scatter(dailyRE['day'],dailyRE['runoff'].replace(0,np.nan),marker='o', c='mediumblue', label='runoff [mm/d]') #we do not want to see zero runoff values
ax2.scatter(dailyRE['day'],dailyRE['erosion'].replace(0,np.nan),marker='o', c='brown', label='erosion [t/ha/d]')
ax2.set_ylabel('daily runoff [mm/day] and erosion [t/ha/day]', fontsize=15)
ax2.tick_params(labelsize=15)
ax2.set_ylim(0,0.5+max(max(dailyRE['runoff']),max(dailyRE['erosion'])))
ax2.legend(loc='upper right', fontsize=15)
fig.tight_layout()
fig.show()

