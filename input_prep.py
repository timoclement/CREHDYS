# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 11:43:11 2023

@author: ticlement
"""
import numpy as np
import os
import csv
import matplotlib.pyplot as plt
import pandas as pd


### GEN ###
def prep_gen(current_wd,plot_name,simdur,dt,nsoil,tp,fc,wp,ksat_sub,df,cl,sa,st,wcf,vfs,rrch,mnch,albed,dx,ninter,nmaize,startdayinter,startdaymaize,dest,bur,wheelinter,wheelmaize,rrini_inter,rrini_maize,kinio,kinim,rrksvarint,rrksvarmaize,kextint,kextmaize,pcoverint,kcovermaize,rbiomcover):
    gen_inp = [
        [plot_name,'','','','Plot name'], #3 spaces between values and explanation so that we have 3 tabs for better reading convenience
        [simdur,'','','','simdur [days]'],
        [dt,'','','','dt [s]'],
        [nsoil,'','','','nsoil'],
        [tp,fc,wp,ksat_sub,df,'','','','tp [m3pores/m3soil]','fc [m3water/m3pores]','tp [m3water/m3pores]','ksat sub [mm/h]','df [cm]'],
        [cl,sa,st,wcf,vfs,'','','','clay [%]','sand','silt','coarse fragments','very fine sand'],
        [rrch,mnch,'','','','RR in channel [mm]', 'Mannings n in channel'],
        [albed,'','','','soil albedo'],
        [dx,'','','','dx [m]'],
        [ninter,nmaize,'','','','number of intercrops periods','number of maize periods'],
        startdayinter + ['','',''] + ['starting day(s) of intercrop periods(s)'],
        [dest,bur,'','','','destruction date','burial date'],
        startdaymaize + ['','','',] + ['starting day(s) of maize period(s)'],
        wheelinter + wheelmaize + ['','',''] + ['number of wheel tracks in intercrop','number of wheel tracks in maize'],
        rrini_inter + ['','',''] + ['initial RR(s) for intercrop period(s)'],
        rrini_maize + ['','',''] + ['initial RR(s) for maize period(s)'],
        [kinio,kinim,'','','','initial Ksat [mm/h] for intercropping and maize periods respectively'],
        rrksvarint + ['','',''] + ['dummy indicating if RR and Ksat are varying (=1) or not (=0) during each intercropping period'],
        rrksvarmaize + ['','',''] + ['dummy indicating if RR and Ksat are varying during each maize period'],
        [kextint,kextmaize,'','','','cover extinction coefficient [m-2/m²] in eq LAI=f(cov) for intercrop and maize respectively'],
        [pcoverint,kcovermaize,'','','','growth rate coefficients for intercrop and maize respectively'],
        [rbiomcover,'','','','coefficient for cover crop cover=f(biom)']
    ]
    
    with open(os.path.join(current_wd,'gen.inp'), 'w', newline='') as csv_file:
        csv_writer=csv.writer(csv_file, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='\\') #delimiter \t means a tab
        csv_writer.writerows(gen_inp)

    print("\033[1;32mgen.inp successfully generated in directory :", current_wd,"\033[0m")


### PARAMHYD ###
def prep_paramhyd(current_wd,kso2,kcho2,kchm2,hfront,c_kso,c_rr,mann,MB_stor,prop_nodep,connectresh):
    paramhyd_inp = [
        [kso2,kcho2,kchm2,hfront,c_kso,c_rr,mann,MB_stor,prop_nodep,connectresh,'','','','final Ksat [mm/h] of overland flow cells','Ksat [mm/h] of wheel tracks in intercropping period','Ksat [mm/h] of wheel tracks in maize period','soil matric potential at wetting front [mm]','soil stability factor [m²/J] of overland flow cells Ksat evolution','soil stability factor of overland flow cell RR evolution','Mannings n for overland flow cells [m^-1/3 / s]','surface storage in microdepressions due to micro-basin tillage [mm]','proportion of surface that is not under microdepression [/]','connectivity threshold [/]']
    ]
    
    with open(os.path.join(current_wd,'paramhyd.inp'), 'w', newline='') as csv_file:
        csv_writer=csv.writer(csv_file, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='\\') #delimiter \t means a tab
        csv_writer.writerows(paramhyd_inp)
        
    print("\033[1;32mparamhyd.inp successfully generated in directory :", current_wd,"\033[0m")


### PARAMEROS ###
def prep_parameros(current_wd,as_,d50,coh,cohw):
    parameros_inp = [
        [as_,d50,coh,cohw,'','','','aggregate stability (Low test) [# drops]','median soil particle diameter [µm]','soil cohesion (torvane) of overland flow cells [kPa]','soil cohesion (torvane) of wheel tracks cells']
    ]
    
    with open(os.path.join(current_wd,'parameros.inp'), 'w', newline='') as csv_file:
        csv_writer=csv.writer(csv_file, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='\\') #delimiter \t means a tab
        csv_writer.writerows(parameros_inp)
        
    print("\033[1;32mparameros.inp successfully generated in directory :", current_wd,"\033[0m")

### PLOT : RECTANGULAR HOMOGENEOUS PLOT ###
def prep_rec_hom_plot(current_wd,dx,xlength,ywidth,slope,yvert_wt,wdth_wt):
    #'''
    if xlength % dx !=0:
        raise ValueError(f"xlength={xlength}m is not a multiple of dx={dx}m")
    if ywidth % dx !=0:
        raise ValueError(f"ywidth={ywidth}m is not a multiple of dx={dx}m")
    if wdth_wt < dx:
        raise ValueError(f"wheel track width={wdth_wt}m is lower than dx={dx}m")
    if wdth_wt % dx !=0:
        raise ValueError(f"wheel track width={wdth_wt}m is not a multiple of dx={dx}m")
    if (np.any(np.array(yvert_wt)%(dx/2))):
        raise ValueError("At least one of vertical wheel track y center position is not located in the middle of a cell or at intersection between two cells.") 
    if (np.any(np.array(yvert_wt)%(dx)==0) and (wdth_wt/dx)%2!=0): #wheel tracks are located between adjacent cells BUT wheel track width is unpair 
        raise ValueError("At least one wheel track is located at intersection between two cells, but the wheel track width does not correspond to an even (pair) number of cells")
    if (np.any(np.array(yvert_wt)%(dx)!=0) and np.any(np.array(yvert_wt)%(dx/2)==0) and (wdth_wt/dx)%2==0): #wheel tracks are located at middle of cell BUT wheel track width is pair 
        raise ValueError("At least one wheel track is located at middle of a cell, but the wheel track width does not correspond to an odd (unpair) number of cells")
    #'''
    flowdir=4 #flow diretion=south=4 in the CREHDYS coding
    nrows=int(xlength/dx) ; ncols=int(ywidth/dx)
    if (nrows*ncols>100485) : raise ValueError(f"The total number of cells {nrows}*{ncols} exceeds the maximum size of spatial vectors in CREHDYS : 100485")
    plot_inp=[[nrows*ncols,nrows,ncols]]
    k=1 #begin at 1 because line 0 is alreay used for ncells, nrows, and ncols
    if yvert_wt: wt=0 #dummy indicating if there are wheel tracks, and on which wheel track of yvert_wt the loop is working on
    else: wt=-1 #no wheel tracks
    
    for i in range(1,nrows+1): #nrows+1 because this means that the loop goes to nrows
        for j in range (1,ncols+1):
            if (wt==-1) : #no wheel tracks at all
                plot_inp.append([k,i,j,slope,flowdir,i,1,1,0]) # 1,1 are soil and crop indice. last 0 indicate that this is NOT a wheel track
                k=k+1
            elif (wt>=0 and ((j*dx<=(yvert_wt[wt]-(wdth_wt/2))) or (j*dx>(yvert_wt[wt]+(wdth_wt/2))))) : # there are wheel tracks but we are not on one (either before or after)
                plot_inp.append([k,i,j,slope,flowdir,i,1,1,0])
                k=k+1
            elif (wt>=0 and j*dx>(yvert_wt[wt]-(wdth_wt/2)) and j*dx<=(yvert_wt[wt]+(wdth_wt/2))) : #there are wheel tracks and we are on one   
                plot_inp.append([k,i,j,slope,flowdir,i,1,1,1]) #1,1 are soil and crop indice. last 1 indicate that this IS a wheel track
                k=k+1
                if (j*dx==(yvert_wt[wt]+(wdth_wt/2)) and wt<len(yvert_wt)-1) : wt=wt+1 # next wheel track
                elif (j*dx==(yvert_wt[wt]+(wdth_wt/2)) and wt==len(yvert_wt)-1) : wt=0 #there are no more wheel tracks on this row after that
                
    with open(os.path.join(current_wd,'plot.inp'), 'w', newline='') as csv_file:
        csv_writer=csv.writer(csv_file, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='\\') #delimiter \t means a tab
        csv_writer.writerows(plot_inp)
        
    print("\033[1;32mplot.inp successfully generated in directory :", current_wd,"\033[0m")
    
    #plotter le plot (field)#
    xcell=[row[1] for row in plot_inp[1:]] #plot_inp[1:] because we do not take line0 which is the ncells nrwos and ncols
    xcoord=[(element * (dx)) - (dx/2) for element in xcell]
    ycell=[row[2] for row in plot_inp[1:]]
    ycoord=[(element * (dx)) - (dx/2) for element in ycell]
    slope=[row[3] for row in plot_inp[1:]]
    flacc=[row[5] for row in plot_inp[1:]]
    channel=[row[8] for row in plot_inp[1:]]
    class_colors={0: 'saddlebrown', 1: 'red', -9999: 'black'}
    colors=[class_colors[value] for value in channel]
    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    ax.scatter(ycoord,xcoord,flacc,c=colors,marker='s',alpha=1)
    ax.set_xlabel('y coord (width) [m]') ; ax.set_ylabel('xcoord (length) [m]') ; ax.set_zlabel('flow acc')
    plt.show
    


### THETA ###
def prep_theta(current_wd,dx,xlength,ywidth,theta_ini):
    nrows=int(xlength/dx) ; ncols=int(ywidth/dx)
    theta_inp=[theta_ini] * (nrows*ncols)
    
    with open(os.path.join(current_wd,'theta.inp'), 'w', newline='') as csv_file:
        csv_writer=csv.writer(csv_file, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='\\') #delimiter \t means a tab
        for value in theta_inp : csv_writer.writerow([value]) #writerow au singulier (pas writerowS!!) car une seule ligne
        
    print("\033[1;32mtheta.inp successfully generated in directory :", current_wd,"\033[0m")



### WEATHER ###
def prep_weather(current_wd,simdur,daily_weather,event_rain):
    daily_weather=pd.read_excel(daily_weather)
    
    event_rain=pd.read_excel(event_rain)
    from datetime import datetime
    
    
    if daily_weather.iloc[:,0].dtype!='datetime64[ns]':
        raise ValueError("type of mdate column in daily weather file is not datetime64[ns]")
    if event_rain.iloc[:,0].dtype!='datetime64[ns]':
        raise ValueError("type of mdate column in daily weather file is not datetime64[ns]")
    if daily_weather.iloc[0,0].date()!=event_rain.iloc[0,0].date():
        raise ValueError(f"time start of daily weather file {daily_weather.iloc[0,0]} is not the same as time start of event rain file {event_rain.iloc[0,0]}")
    if daily_weather.iloc[-1,0].date()!=event_rain.iloc[-2,0].date():
        raise ValueError(f"time end of daily weather file {daily_weather.iloc[-1,0]} is not the same as time end of event rain file {event_rain.iloc[-2,0]}")
    if int((daily_weather.iloc[-1,0].date() - daily_weather.iloc[0,0].date()).days)+1 != simdur:
        raise ValueError(f"time duration of input weather files (daily and event) does not correspond to simdur={simdur} days")

    ievent=0 ; jevent=0 #indices for event input file reading
    weather_inp = []
    
    for day in range (1,simdur+1) :
        pluday=0
        jevent=ievent
        starteventk=-9999
        while event_rain.iloc[ievent,0].date() == daily_weather.iloc[day-1,0].date() : #read total rainfall of this day
            if (starteventk==-9999 and event_rain.iloc[ievent,1]>0) : starteventk = ievent # indice of event start
            pluday=pluday+event_rain.iloc[ievent,1]
            ievent=ievent+1 #At end of the while loop, ievent=first indice of the next day in event input file
        
        if (pluday==0) : # NO RAIN THIS DAY  
             weather_inp.append([day, daily_weather.iloc[day-1,3],daily_weather.iloc[day-1,6],round(daily_weather.iloc[day-1,10],2),0,0]) # write daily weather, with dummy=0 indicating that there is no rain event this day 

        if (pluday>0) : # RAIN THIS DAY
            starteventhour=round((event_rain.iloc[starteventk,0] - event_rain.iloc[starteventk,0].replace(hour=0,minute=0,second=0)).total_seconds()/3600 , 2)
            weather_inp.append([day, daily_weather.iloc[day-1,3],daily_weather.iloc[day-1,6],round(daily_weather.iloc[day-1,10],2),1,starteventhour])# write daily weather, with dummy=1 indicating that there is a rain event this day
            pluday2=0 #daily cumulative rainfall to track end of event
            tevent=0
            weather_inp.append([day,0,0,0,0,0]) # write first zero event rain
            while (event_rain.iloc[jevent,0].date() == daily_weather.iloc[day-1,0].date() ) : # while in same day
                deltatmin = round((event_rain.iloc[jevent,0] - event_rain.iloc[jevent-1,0]).total_seconds()/60 , 0) # delta time between two rain records [min]
                rainrate =  round(event_rain.iloc[jevent,1] / (deltatmin/60) , 1) #rainfall rate [mm/h] during these two rain records
                if (jevent>=starteventk and pluday2<pluday):
                    tevent = tevent + deltatmin
                    weather_inp.append([day, tevent, rainrate, 0,0,0]) # write event rain
                #if (tevent==endeventtime*60) : weather_inp.append([day, int(endeventtime*60-starteventtime*60+deltatmin), 0, 0,0,0]) #add last record with zero rain
                pluday2=pluday2+event_rain.iloc[jevent,1]
                jevent=jevent+1
            
            weather_inp[len(weather_inp)-1][4]=1 #last record of the event -> dummy=1 indicating that it is end of rain event
                
    
    
    with open(os.path.join(current_wd,'weather.inp'), 'w', newline='') as csv_file:
            csv_writer=csv.writer(csv_file, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='\\') #delimiter \t means a tab
            csv_writer.writerows(weather_inp)
    
    print("\033[1;32mweather.inp successfully generated in directory :", current_wd,"\033[0m")

        


