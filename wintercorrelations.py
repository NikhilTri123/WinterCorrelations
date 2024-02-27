import csv
import numpy as np
import scipy.stats as stats
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import warnings
from shapely.errors import ShapelyDeprecationWarning

warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)


def getCorrelation(data, snowVals, month):
    dim = data.sel(time=np.datetime64('2017-01-01'))
    all_heights = []
    used_snow_vals = []
    for i in range(len(snowVals)):
        if snowVals[i] == 'T':
            snowVals[i] = 0
        try:
            value = float(snowVals[i])
            used_snow_vals.append(value)
            height_year = data.sel(time=f"{1960 + i}-{format(month, '02d')}-01").values.flatten()
            all_heights.append(height_year)
        except ValueError:
            continue
    all_heights = np.array(all_heights)

    flipped_heights = []
    for j in range(len(all_heights[0])):
        flipped_heights.append(all_heights[:, j])
    flipped_heights = np.array(flipped_heights)

    corrList = []
    sigList = []
    for k in range(len(flipped_heights)):
        corr = stats.pearsonr(flipped_heights[k], used_snow_vals)
        corrList.append(corr[0])
        sigList.append(corr[1])
    print("Correlation created")
    corrList = np.reshape(corrList, (len(dim), len(dim[0])))
    sigList = np.reshape(sigList, (len(dim), len(dim[0])))
    return corrList, sigList


snow_data = "iad_snowfall.csv"
snow_list = list(csv.reader(open(snow_data)))
months = snow_list[0]
snow_data = np.array(snow_list[1:])
print(snow_list)

hgt_file = "hgtwintercoarse.nc"
hgt_data = xr.open_dataset(hgt_file).z
print(hgt_data)

mon_data = snow_data[:, months.index(" Feb")]
print(mon_data)
corr_data, sig_data = getCorrelation(hgt_data, mon_data, 2)

corr_dataset = xr.Dataset(
    data_vars=dict(
        decComp=(["latitude", "longitude"], corr_data)
    ),
    coords=dict(
        latitude=(["latitude"], hgt_data.latitude.values),
        longitude=(["longitude"], hgt_data.longitude.values),
    )
)
corr_data_array = corr_dataset.decComp.sel(latitude=slice(90, -5), longitude=slice(160, 360))
print(corr_data_array)

sig_dataset = xr.Dataset(
    data_vars=dict(
        decComp=(["latitude", "longitude"], sig_data)
    ),
    coords=dict(
        latitude=(["latitude"], hgt_data.latitude.values),
        longitude=(["longitude"], hgt_data.longitude.values),
    )
)
sig_data_array = sig_dataset.decComp.sel(latitude=slice(90, -5), longitude=slice(160, 360))

# plot cartopy map and various features
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111, projection=ccrs.LambertConformal(central_longitude=-100))
ax.add_feature(cf.LAND)
ax.add_feature(cf.STATES, linewidth=0.2, edgecolor="gray")
ax.add_feature(cf.BORDERS, linewidth=0.3)
ax.coastlines(linewidth=0.5, resolution='50m')
ax.set_extent([-145, -55, 7, 74])

# plot gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), draw_labels=True, linewidth=1, color='gray', alpha=0.5,
                  linestyle='--')
gl.top_labels = gl.right_labels = False
gl.xlabel_style = {'size': 6, 'weight': 'bold', 'color': 'gray'}
gl.ylabel_style = {'size': 6, 'weight': 'bold', 'color': 'gray'}

# add data and colorbar
plt.contourf(corr_data_array.longitude, corr_data_array.latitude, corr_data_array, 60, levels=np.arange(-0.5, 0.5, 0.01),
             extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
cbar = plt.colorbar(pad=0.017, aspect=25, shrink=1)
cbar.ax.tick_params(labelsize=8)

# add significance contouring
contour = plt.contour(sig_data_array.longitude, sig_data_array.latitude, sig_data_array, 60, levels=[0, 0.05],
                      transform=ccrs.PlateCarree(), colors='black', linestyles='dashed', linewidths=0)
for col in contour.collections:
    col.set_hatch('xx')
    col.set_alpha(0.2)

# add titling
mainTitle = "IAD Feb Snow Correlated to 500mb Heights"
subTitle = "\nData from ERA5 and xmACIS2"
plt.title(mainTitle + subTitle, fontsize=9, weight='bold', loc='left')
plt.title("DCAreaWx", fontsize=9, weight='bold', loc='right', color='gray')
ax.text(-141.5, 1.5, 'Hatched regions indicate significant correlations', fontsize=9, weight='bold',
        transform=ccrs.PlateCarree())

plt.savefig(r"jan_snow_corr.png", dpi=300, bbox_inches='tight')
plt.show()

# ----------------------------------------------------------------------------------------------------------------------

# # plot cartopy map and various features
# plt.figure(figsize=(12, 6))
# ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
# ax.add_feature(cf.LAND)
# ax.add_feature(cf.STATES, linewidth=0.2, edgecolor="gray")
# ax.add_feature(cf.BORDERS, linewidth=0.3)
# ax.coastlines(linewidth=0.5, resolution='50m')
#
# # plot gridlines
# gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), draw_labels=True, linewidth=1, color='gray', alpha=0.5,
#                   linestyle='--')
# gl.top_labels = gl.right_labels = False
# gl.xlabel_style = {'size': 6, 'weight': 'bold', 'color': 'gray'}
# gl.ylabel_style = {'size': 6, 'weight': 'bold', 'color': 'gray'}
#
# # add data and colorbar
# plt.contourf(corr_data_array.longitude, corr_data_array.latitude, corr_data_array, 60, levels=np.arange(-0.5, 0.5, 0.01),
#              extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
# cbar = plt.colorbar(pad=0.017, aspect=25, shrink=0.78)
# cbar.ax.tick_params(labelsize=6)
#
# # add significance contouring
# contour = plt.contour(sig_data_array.longitude, sig_data_array.latitude, sig_data_array, 60, levels=[0, 0.05],
#                       transform=ccrs.PlateCarree(), colors='black', linestyles='dashed', linewidths=0)
# for col in contour.collections:
#     col.set_hatch('xx')
#     col.set_alpha(0.2)
#
# # add titling
# mainTitle = "IAD Feb Snow Correlated to 500mb Heights"
# subTitle = "\nData from ERA5 and xmACIS2"
# plt.title(mainTitle + subTitle, fontsize=9, weight='bold', loc='left')
# plt.title("DCAreaWx", fontsize=9, weight='bold', loc='right', color='gray')
# ax.text(172, 2, 'Hatched regions indicate significant correlations', fontsize=9, weight='bold',
#         transform=ccrs.PlateCarree())
#
# plt.savefig(r"jan_snow_corr.png", dpi=300, bbox_inches='tight')
# plt.show()
