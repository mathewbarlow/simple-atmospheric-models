# coded in Python 3.7
# but in a very Fortran style (may cause nausea and uncontrollable crying)
#
# Simple program to solve trajectories for inertial oscillations
#
# This version makes an animation, in gif format
#
# Note that interpretation of physics is complicated; e.g.,
# Durran, D.R., 1993: Is the Coriolis Force Really Responsible for the Inertial 
# Oscillation?. Bull. Amer. Meteor. Soc., 74, 2179–2184, 
# https://doi.org/10.1175/1520-0477(1993)074<2179:ITCFRR>2.0.CO;2 
#
# Coded by M. Barlow, please email questions or comments to 
# Mathew_Barlow@uml.edu
#
# initially released 24 Oct 2019
# last modified 15 Nov 2019 (analytic sol'n added, animation added)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import imageio
import os

# user-defined values
nt = 60*24*3    # number of time steps
dt = 60    # time step, in seconds
lon0 = 288.9   # starting longitude, in degrees
lat0 = 42.4    # starting latitude, in degrees
u0 = 20       # initial zonal wind, in m/s
v0 = 0         # initial meridional wind, in m/s

# constants
d2r = np.pi/180   # conversion factor for degrees to radians
omega = 7.292E-5  # angular velocity of earth
re = 6.371E6      # radius of Earth, in meters

# these are needed for first time step
xold = 0.0    # initial x position, in meters, relative to starting longitude
yold = 0.0    # initial y position, in meters, relative to starting latitude
uold = u0     # initial zonal wind
vold = v0     # initial meridional wind
lonold = lon0
latold = lat0

# arrays to store trajectory information at each time step
x_all = np.zeros(nt+1)
y_all = np.zeros(nt+1)
lon_all = np.zeros(nt+1)
lat_all = np.zeros(nt+1)
u_all = np.zeros(nt+1)
v_all = np.zeros(nt+1)
pairs = np.zeros((nt+1,2))

# store initial conditions
x_all[0] = xold
y_all[0] = yold
lon_all[0] = lonold
lat_all[0] = latold
u_all[0] = uold
v_all[0] = vold

# analytic solution for fixed f and dx
t_all = np.arange(0,nt*dt,dt)
f0 = 2*omega*np.sin(lat0*d2r)
x_lytic_all = u0*np.sin(f0*t_all)/f0 - v0*np.cos(f0*t_all)/f0 + v0/f0
y_lytic_all = v0*np.sin(f0*t_all)/f0 + u0*np.cos(f0*t_all)/f0 - u0/f0
lon_lytic_all = lon0 + x_lytic_all*np.cos(lat0*d2r)*360/(2*np.pi*re)
lat_lytic_all = lat0 + y_lytic_all*360/(2*np.pi*re)
u_lytic_all = u0*np.cos(f0*t_all) + v0*np.sin(f0*t_all)
v_lytic_all = v0*np.cos(f0*t_all) - u0*np.sin(f0*t_all)


# Note: always keeping track of three times:  old, current, and new. 
# For example:  xold, x, and xnew.

#First time step is forward difference
x = xold + dt*uold
y = yold + 2*dt*vold
f = 2*omega*np.sin(latold*d2r)
u = uold + dt*f*vold
v = vold - dt*f*uold

lat = latold + (y-yold)*360/(2*np.pi*re)
lon = lonold + (x-xold)*np.cos(lat*d2r)*360/(2*np.pi*re)

#   store values 
x_all[1] = x
y_all[1] = y
lon_all[1] = lon
lat_all[1] = lat
u_all[1] = u
v_all[1] = v


#All subsequent steps are center difference
it = 2
while it<=nt:

#   calculate new values    
    xnew = xold + 2*dt*u
    ynew = yold + 2*dt*v
    f = 2*omega*np.sin(lat*d2r)
    unew = uold + 2*dt*f*v
    vnew = vold - 2*dt*f*u
    
    lonnew = lonold + (xnew-xold)*np.cos(lat*d2r)*360/(2*np.pi*re)
    latnew = latold + (ynew-yold)*360/(2*np.pi*re)

#   update values for next time step    
    xold = x
    yold = y
    x = xnew
    y = ynew
    uold = u
    vold = v
    u = unew
    v = vnew
    lonold = lon
    latold = lat
    lon = lonnew
    lat = latnew
    
#   store values 
    x_all[it] = x
    y_all[it] = y
    lon_all[it] = lon
    lat_all[it] = lat
    u_all[it] = u
    v_all[it] = v
    
#    print('time step: ',it, 'x: ', x, 'y: ', y)
    
    it = it + 1
    
    
# plot
plt.clf()

plt.figure(1)

fudge=1E-2
lon1 = np.int(lon0)-5
lon2 = np.int(lon0)+5
lat1 = np.int(lat0)-5
lat2 = np.int(lat0)+5

ax = plt.axes(projection=ccrs.PlateCarree( ))
ax.set_extent([lon1-fudge,lon2+fudge,lat1-fudge,lat2+fudge],crs=ccrs.PlateCarree())

gl = ax.gridlines(draw_labels=True)
ax.coastlines('50m', linewidth=0.8,color='gray')
states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                             name='admin_1_states_provinces_shp')
ax.add_feature(states, linewidth=0.5,edgecolor='grey')

plt.plot(lon_all,lat_all,'-',markersize=2,transform=ccrs.PlateCarree())
plt.plot(lon_all[-1],lat_all[-1],'o',markersize=10,color='darkblue',
transform=ccrs.PlateCarree())

gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
dlon=(lon2-lon1)/5
dlat=(lat2-lat1)/5
gl.xlocator = mticker.FixedLocator(np.arange(lon1-360,lon2+dlon-360,dlon))
gl.ylocator = mticker.FixedLocator(np.arange(lat1,lat2+dlat,dlat))
plt.title('Trajectory')


plt.figure(2)
plt.ylim([199.75, 200.25])
plt.plot((u_all*u_all+v_all*v_all)/2)
plt.title('Kinetic Energy')


plt.figure(3)

fudge=1E-2
lon1 = np.int(lon0)-5
lon2 = np.int(lon0)+5
lat1 = np.int(lat0)-5
lat2 = np.int(lat0)+5

ax = plt.axes(projection=ccrs.PlateCarree( ))
ax.set_extent([lon1-fudge,lon2+fudge,lat1-fudge,lat2+fudge],crs=ccrs.PlateCarree())

gl = ax.gridlines(draw_labels=True)
ax.coastlines('50m', linewidth=0.8,color='gray')
states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                             name='admin_1_states_provinces_shp')
ax.add_feature(states, linewidth=0.5,edgecolor='grey')

plt.plot(lon_lytic_all,lat_lytic_all,'-',markersize=2,color='red',
transform=ccrs.PlateCarree())
plt.plot(lon_all,lat_all,'-',markersize=2,transform=ccrs.PlateCarree())
plt.plot(lon_all[-1],lat_all[-1],'o',markersize=10,color='darkblue',
transform=ccrs.PlateCarree())

gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
dlon=(lon2-lon1)/5
dlat=(lat2-lat1)/5
gl.xlocator = mticker.FixedLocator(np.arange(lon1-360,lon2+dlon-360,dlon))
gl.ylocator = mticker.FixedLocator(np.arange(lat1,lat2+dlat,dlat))
plt.title('Trajectory')


plt.close(fig='all')

fudge=1E-2
lon1 = np.int(lon0)-5
lon2 = np.int(lon0)+5
lat1 = np.int(lat0)-5
lat2 = np.int(lat0)+5
states = NaturalEarthFeature(category='cultural', scale='50m', 
         facecolor='none', name='admin_1_states_provinces_shp')
dlon=(lon2-lon1)/5
dlat=(lat2-lat1)/5

images = []
filenames = []

nframe=90
iframe=0
while iframe<=nframe:
    
    it=np.int(iframe*(nt-1)/nframe)+1
    plt.figure(iframe,figsize=plt.figaspect(0.5))
    
    ax = plt.axes(projection=ccrs.PlateCarree( ))
    ax.set_extent([lon1-fudge,lon2+fudge,lat1-fudge,
                  lat2+fudge],crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True)
    ax.coastlines('50m', linewidth=0.8,color='gray')
    ax.add_feature(states, linewidth=0.5,edgecolor='grey')
    
    plt.plot(lon_all[0:it],lat_all[0:it],'-',markersize=2,
             transform=ccrs.PlateCarree())
    plt.plot(lon_all[it],lat_all[it],'o',markersize=10,color='darkblue',
             transform=ccrs.PlateCarree())
    
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(lon1-360,lon2+dlon-360,dlon))
    gl.ylocator = mticker.FixedLocator(np.arange(lat1,lat2+dlat,dlat))
    plt.title('Trajectory')
    
    mydir='/Users/mathew_barlow/downloads/'
    filename=mydir +'temp'+ '{:04d}'.format(iframe)+'.png'
    plt.savefig(filename,bbox_inches='tight', dpi=150)
    
    images.append(imageio.imread(filename))
    filenames.append(filename)
                    
    iframe=iframe+1

plt.clf()

imageio.mimsave(mydir+'inert.gif', images)
for file in filenames:
    os.remove(file)
