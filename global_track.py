# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 15:40:33 2017

@author: bling
"""

import sys
#import pytz
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from global_track_functions import get_drifter,get_fvcom,draw_basemap
from matplotlib import animation

 
track_days = 3
MODEL = 'global' #'global'GOM3
track_direction = 'forword' # backword,forward

points = {'lats':[],'lons':[]}  # collect boundary points.
drifter_points = dict(lon=[],lat=[])
model_points = dict(lon=[],lat=[])

######################################## Drifter ##################################
print 'Drifter parts.'
drifter_ID = '160420692'#152300811 
# if raw data, use "drift_X.dat";if want to get drifter data in database, use "None"
INPUT_DATA = 'drift_X.dat'#'drift_jml_2015_1.dat'
start_time = datetime.utcnow()-timedelta(track_days) #datetime(2015,1,24,0,0,0,0,pytz.UTC)
drifter = get_drifter(drifter_ID, INPUT_DATA)
dr_points = drifter.get_track(start_time,track_days)
drifter_points['lon'].extend(dr_points['lon']); drifter_points['lat'].extend(dr_points['lat'])
print "drifter points: ",len(dr_points['lon']),'\nlast point(',dr_points['lat'][-1],',',dr_points['lon'][-1],')'
#np.savez('drifter_points.npz',lon=drifter_points['lon'],lat=drifter_points['lat'])#'''

####################################### Model ###################################
print 'Model parts.'
centerpoint = (dr_points['lat'][-1],dr_points['lon'][-1])
#centerpoint = (7.77527332,52.91501617)
bordersidele = 2
start_time = dr_points['time'][-1]-timedelta(track_days)
#start_time = datetime.utcnow()
end_time = start_time + timedelta(track_days)
print 'Start time: ',start_time,'\n','End time: ',end_time 

get_obj =  get_fvcom(MODEL)
modeltime = get_obj.get_url(start_time,end_time)
psqus = get_obj.get_data(centerpoint,bordersidele) # b_points is model boundary points.
points['lons'].extend(psqus[1]);points['lats'].extend(psqus[0])
point = get_obj.get_track(dr_points['lon'][0],dr_points['lat'][0],track_direction)
#point = get_obj.get_track(52.91501617,7.77527332,track_direction)
model_points['lon'].extend(point['lon']); model_points['lat'].extend(point['lat'])
print point
###################################### Plot ########################################
fig = plt.figure() #figsize=(16,9)
ax = fig.add_subplot(111)
#draw_basemap(ax, points)  # points is using here
ax.plot(drifter_points['lon'],drifter_points['lat'],'b')
ax.plot(model_points['lon'],model_points['lat'],'r')
#ax.plot()
plt.show()
