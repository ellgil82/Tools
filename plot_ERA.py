import iris
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import cartopy
import os
import cartopy.crs as ccrs
import sys
sys.path.append('/users/ellgil82/scripts/Tools/')
from divg_temp_colourmap import shiftedColorMap

def load_files(case): #period should be either 'pre' for pre-foehn, or 'onset' for foehn conditions
    if case == 'CS1':
        os.chdir('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS1/')  # path to data
        surf = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS1/20160510T1200Z_sfc_temp_MSLP.nc'
        winds = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS1/20160510T1200Z_750_wind_geopot.nc'
    elif case == 'CS2':
        os.chdir('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/')
        surf = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/20160526T0000Z_sfc_temp_MSLP.nc'
        winds = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/20160526T0000Z_750_wind_geopot.nc'
    print ('importing cubes...')
    u = iris.load_cube(winds, 'eastward_wind')
    v = iris.load_cube(winds, 'northward_wind')
    P = iris.load_cube(surf, 'Mean sea level pressure')
    P.convert_units('hPa')
    geopot = iris.load_cube(winds, 'geopotential')
    air_T = iris.load_cube(surf, '2 metre temperature')
    air_T.convert_units('celsius')
    lon = P.coord('longitude')
    lat = P.coord('latitude')
    return u[0,:,:], v[0,:,:], air_T[0,:,:], geopot[0,:,:], P[0,:,:], lat, lon


def synoptic_means(): # synoptic monthly means at 12Z (t=0) for (Jan) 2011
    synop = '/data/mac/ellgil82/cloud_data/ERA_Int_wind_cloud_monthly_synop_means_at_12Z_2011.nc'
    cubes = iris.load(synop)
    v_syn = cubes [2]
    u_syn = cubes[1]
    cl_syn = cubes[0]
    u_syn = u_syn[0,200:220,390:410]#,1220:1290,2320:2431]
    v_syn = v_syn[0,200:220,390:410]#,1220:1290,2320:2431]
    cl_syn = cl_syn[0,200:220,390:410]#,1220:1290,2320:2431]
    lon = u_syn.coord('longitude')
    lat = u_syn.coord('latitude')
    x,y = np.meshgrid(lon.points, lat.points)
    clevs = np.arange(0,1.01, 0.1)
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())# maximise window
    c = ax.contourf(lon.points, lat.points, cl_syn.data, levels= clevs, cmap='Greys_r', vmin=0.0, vmax=1.0)
    CBAxes = fig.add_axes([0.75, 0.1, 0.03, 0.7])
    cb = plt.colorbar(c, cax = CBAxes)
    cb.ax.tick_params(labelsize=24)
    cb.set_label('Cloud fraction', fontsize = 24, labelpad = 10)
    ax.coastlines(resolution='50m', linewidth=2)
    plt.subplots_adjust(left = 0.05, right = 0.8, top = 0.8, bottom = 0.1)
    q = ax.quiver(x[::3,::3],y[::3,::3],u_syn.data[::3,::3],v_syn.data[::3,::3], pivot='middle',  scale = 100)
    plt.quiverkey(q,0.25, 0.85, 10, r'$10$ $m$ $s^{-1}$', labelpos='N', fontproperties={'size':'24', 'weight':'bold'},
                       coordinates='figure')
    plt.draw()
    plt.rcParams['svg.fonttype'] = 'none'
    plt.savefig('/users/ellgil82/figures/Cloud data/f152/ERA_Int_synoptic_winds_cloud_frac_750hPa_2011_Jan_12Z.eps', dpi = 300)
    plt.savefig('/users/ellgil82/figures/Cloud data/f152/ERA_Int_synoptic_winds_cloud_frac_750hPa_2011_Jan_12Z.png')

## Caption: Synoptic conditions on <date> at 12:00 UTC over the Antarctic Peninsula, derived from high resolution ERA-5 HRES operational analysis
## at 31 km resolution. Colour contours show the 2 m temperature in degrees celsius, while the overlaid vectors show 750 hPa winds.

def case():
    fig = plt.figure(frameon=False, figsize=(10,12)) # !!change figure dimensions when you have a larger model domain
    fig.patch.set_visible(False)
    label_dict = {0:'a', 1:'b'}
    dates = ['20160508T0000Z', '20160526T0000Z']#
    for i in [0]:#,1]:
        ax = fig.add_subplot(len(dates),1,i+1,projection=ccrs.PlateCarree())
        surf_file = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/'+dates[i]+'_sfc_temp_MSLP.nc'
        wind_file = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/'+dates[i]+'_750_wind_geopot.nc'
        #synop_file = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/'+dates[i]+'_synop_surf_era5.nc'
        u = iris.load_cube(wind_file, 'eastward_wind')
        v = iris.load_cube(wind_file, 'northward_wind')
        P = iris.load_cube(surf_file, 'air_pressure_at_sea_level')
        P.convert_units('hPa')
        geopot = iris.load_cube(wind_file, 'geopotential')
        geopot = geopot[0,:,:]
        P = P[0,:,:]
        u = u[0,:,:]
        v = v[0,:,:]
        sfc_T = iris.load_cube(synop_file, '2 metre temperature')
        sfc_T.convert_units('celsius')
        lon = P.coord('longitude')
        lat = P.coord('latitude')
        x,y = np.meshgrid(lon.points, lat.points)
        clevs = np.arange(np.floor(np.min(sfc_T.data))-3, np.ceil(np.max(sfc_T.data))+3)
        ax.set_extent((min(lon.points), max(lon.points), min(lat.points), max(lat.points)))
        ax.axis('off')
        ax.tick_params(axis='both', which='both', length=0, labelbottom='off', labelleft='off')
        CMap = shiftedColorMap(cmap=matplotlib.cm.bwr, min_val=np.floor(np.min(sfc_T.data)), max_val=np.ceil(np.max(sfc_T.data)), name='CMap', var = sfc_T)
        ax.coastlines(resolution='50m', linewidth=2)
        c = ax.contourf(lon.points, lat.points, sfc_T[0,:,:].data, cmap=CMap, levels=clevs, transform=ccrs.PlateCarree(), latlon = True)#  vmin=np.floor(np.min(sfc_T.data)), vmax=np.ceil(np.max(sfc_T.data)),
        geo = ax.contour(lon.points, lat.points, geopot.data, colors = '0.3',  linewidths=3, latlon = True)
        MSLP = ax.contour(lon.points, lat.points, P.data, colors = '0.3',  linewidths=3, levels = range(996,1016,4) , latlon = True)
        ax.clabel(MSLP, v=[960, 968, 976, 982, 990], inline=True, inline_spacing=3, fontsize=28, fmt='%1.0f')
        q = ax.quiver(x[::20, ::20], y[::20, ::20], u.data[::20, ::20], v.data[::20, ::20],
                         pivot='middle')  # ,  scale = 200)
        ax.text(x= -73, y=-62, s=label_dict[i], fontsize=32, fontweight='bold', color='k')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
    CBAxes = fig.add_axes([0.25, 0.12, 0.5, 0.03])
    cb = plt.colorbar(c, cax=CBAxes, orientation = 'horizontal')
    cb.ax.tick_params(labelsize=24)
    cb.set_label('2 m air temperature ($^\circ$C)', fontsize=24, labelpad=10, color = 'dimgrey')
    cb.outline.set_edgecolor('dimgrey')
    cb.outline.set_linewidth(2)
    cb.ax.tick_params(labelsize=24, labelcolor='dimgrey', pad=10)
    plt.quiverkey(q, 0.2, 0.88, 50, r'$50$ $m$ $s^{-1}$', labelcolor = 'dimgrey', color = 'dimgrey',  labelpos='N', fontproperties={'size': '24', 'weight': 'bold'},
                  coordinates='figure')
    plt.draw()
    plt.rcParams['svg.fonttype'] = 'none'
    [l.set_visible(False) for (w, l) in enumerate(cb.ax.xaxis.get_ticklabels()) if w % 2 != 0]
    plt.subplots_adjust(left = 0.1, bottom = 0.2, right = 0.9, top = 0.85)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/ERA5_winds_750hPa_2m_temp_CS1_CS2.eps')
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/ERA5_winds_750hPa_2m_temp_CS1_CS2.png')
    plt.show()


#case()
#synoptic_means()

def both():
    fig, axs = plt.subplots(1, 2, figsize=(20, 12), frameon=False, subplot_kw={'projection': ccrs.PlateCarree()})
    axs.flatten()
    lab_dict = {0: 'a', 1: 'b'}
    # ax = fig.add_axes([0.18, 0.25, 0.75, 0.63], frameon=False) # , projection=ccrs.PlateCarree())#
    case_list = ['CS1', 'CS2']
    for a in [0,1]:
        u_wind, v_wind, air_T, geopot, P, lat, lon = load_files(case = case_list[a])
        x,y = np.meshgrid(lon.points, lat.points)
        bwr_zero = shiftedColorMap(cmap=matplotlib.cm.bwr, min_val=-20., max_val=10., name='bwr_zero', var=air_T.data,
                                   start=0., stop=1.)
        c = axs[a].contourf(lon.points, lat.points, air_T.data, cmap=bwr_zero, vmin=-20., vmax=10.)#, levels=[])
        MSLP = axs[a].contour(lon.points, lat.points, P.data, levels=range(960, 1020, 4), colors='#222222', linewidths=3,)
        geoP = axs[a].contour(lon.points, lat.points, geopot.data, colors='#222222', linewidths=3,)
        axs[a].coastlines(resolution='50m', linewidth=2, color='#535454')
        q = axs[a].quiver(x[::25,::25],y[::25,::25],u_wind.data[::25,::25],v_wind.data[::25,::25], pivot='middle', color='#414345', scale = 200)
        axs[a].text(0.1, 0.85, transform = axs[a].transAxes, s=lab_dict[a], fontsize=32, fontweight='bold', color='dimgrey')
        axs[a].tick_params(which='both', axis='both', labelsize=34, labelcolor='dimgrey', pad=10, size=0, tick1On=False,
                           tick2On=False)
        PlotLonMin = np.min(lon.points)
        PlotLonMax = np.max(lon.points)
        PlotLatMin = np.min(lat.points)
        PlotLatMax = np.max(lat.points)
        XTicks = np.linspace(PlotLonMin, PlotLonMax, 3)
        XTickLabels = [None] * len(XTicks)
        for i, XTick in enumerate(XTicks):
            if XTick < 0:
                XTickLabels[i] = '{:.0f}{:s}'.format(np.abs(XTick), '$^{\circ}$W')
            else:
                XTickLabels[i] = '{:.0f}{:s}'.format(np.abs(XTick), '$^{\circ}$E')
        plt.sca(axs[a])
        plt.xticks(XTicks, XTickLabels)
        axs[a].set_xlim(PlotLonMin, PlotLonMax)
        axs[a].tick_params(which='both', pad=10, labelsize=34, color='dimgrey')
        YTicks = np.linspace(PlotLatMin, PlotLatMax, 4)
        YTickLabels = [None] * len(YTicks)
        for i, YTick in enumerate(YTicks):
            if YTick < 0:
                YTickLabels[i] = '{:.0f}{:s}'.format(np.abs(YTick), '$^{\circ}$S')
            else:
                YTickLabels[i] = '{:.0f}{:s}'.format(np.abs(YTick), '$^{\circ}$N')
        plt.sca(axs[a])
        plt.yticks(YTicks, YTickLabels)
        axs[a].set_ylim(PlotLatMin, PlotLatMax)
        axs[a].set_title(case_list[a], fontsize=34, color='dimgrey')
        axs[a].spines['right'].set_visible(False)
        axs[a].spines['left'].set_visible(False)
        axs[a].spines['top'].set_visible(False)
        axs[a].spines['bottom'].set_visible(False)
    CBarXTicks = [-30, -10, 10]  # CLevs[np.arange(0,len(CLevs),int(np.ceil(len(CLevs)/5.)))]
    CBAxes = fig.add_axes([0.35, 0.22, 0.3, 0.03])
    CBar = plt.colorbar(c, cax=CBAxes, orientation='horizontal', ticks=CBarXTicks)  #
    CBar.set_label('2 m air temperature ($^{\circ}$C)', fontsize=34, labelpad=10, color='dimgrey')
    CBar.solids.set_edgecolor("face")
    CBar.outline.set_edgecolor('dimgrey')
    CBar.ax.tick_params(which='both', axis='both', labelsize=34, labelcolor='dimgrey', pad=10, size=0, tick1On=False,
                        tick2On=False)
    CBar.outline.set_linewidth(2)
    plt.sca(axs[1])
    plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
    # yaxis.set_label_coords(1.27, 0.5)
    plt.quiverkey(q, 0.51, 0.9, 20, r'$20$ $m$ $s^{-1}$', labelpos='N', color='#414345', labelcolor='#414345',
                  fontproperties={'size': '32', 'weight': 'bold'},
                  coordinates='figure', )
    plt.draw()
    plt.subplots_adjust(bottom=0.25, top=0.85, left = 0.15, right = 0.85)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/ERA_5_winds_air_T_geopotential.eps')
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/ERA_5_winds_air_T_geopotential.png')
    plt.show()

both()
