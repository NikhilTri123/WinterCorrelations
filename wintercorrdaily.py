import csv
import numpy as np
import scipy.stats as stats
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf


def getCorrelation(data, snowVals, lagtime):
    all_heights = []
    used_snow_vals = []
    for snowVal in snowVals:
        if int(snowVal[0][5:7]) not in [1, 2] or snowVal[1] == ' M':
            continue

        if snowVal[1] == ' T':
            snowVal[1] = 0
        used_snow_vals.append(float(snowVal[1]))
        height_year = data.sel(time=snowVal[0]).values.flatten()
        all_heights.append(height_year)

    all_heights = np.array(all_heights)

    if lagtime != 0:
        all_heights = all_heights[:lagtime]
        used_snow_vals = used_snow_vals[lagtime * -1:]

    corrList = []
    sigList = []
    flipped_heights = all_heights.T
    for pixel in flipped_heights:
        corr = stats.pearsonr(pixel, used_snow_vals)
        corrList.append(corr[0])
        sigList.append(corr[1])
    corrList = np.reshape(corrList, data.shape[1:])
    sigList = np.reshape(sigList, data.shape[1:])
    print("Correlation created")
    return corrList, sigList


def getSnowData(snow_path):
    snow_list = list(csv.reader(open(snow_path)))
    snow_data = np.array(snow_list[1:])
    print(snow_data.shape)
    return snow_data


hgt_file = "hgtwinterdailycoarse.nc"
hgt_data = xr.open_dataset(hgt_file).z.sel(time=slice("1980-01-01", "2022-12-31"), expver=1)
print(hgt_data)


def plot(lagtime, snow_paths):
    print(lagtime)
    all_corr_data = []
    all_sig_data = []
    for snow_path in snow_paths:
        snow_data = getSnowData(snow_path)
        corr_data, sig_data = getCorrelation(hgt_data, snow_data, lagtime)
        all_corr_data.append(corr_data)
        all_sig_data.append(sig_data)
    corr_data = np.mean(all_corr_data, axis=0)
    sig_data = np.mean(all_sig_data, axis=0)

    corr_dataset = xr.Dataset(
        data_vars=dict(
            correlation=(["latitude", "longitude"], corr_data)
        ),
        coords=dict(
            latitude=(["latitude"], hgt_data.latitude.values),
            longitude=(["longitude"], hgt_data.longitude.values),
        )
    )
    corr_data_array = corr_dataset.correlation

    sig_dataset = xr.Dataset(
        data_vars=dict(
            significance=(["latitude", "longitude"], sig_data)
        ),
        coords=dict(
            latitude=(["latitude"], hgt_data.latitude.values),
            longitude=(["longitude"], hgt_data.longitude.values),
        )
    )
    sig_data_array = sig_dataset.significance

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
    plt.contourf(corr_data_array.longitude, corr_data_array.latitude, corr_data_array, 60, levels=np.arange(-0.15, 0.15, 0.005),
                 extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
    cbar = plt.colorbar(pad=0.017, aspect=25, shrink=0.84)
    cbar.ax.tick_params(labelsize=6)

    # add significance contouring
    contour = plt.contour(sig_data_array.longitude, sig_data_array.latitude, sig_data_array, 60, levels=[0, 0.05],
                          transform=ccrs.PlateCarree(), colors='black', linestyles='dashed', linewidths=0)
    for col in contour.collections:
        col.set_hatch('xx')
        col.set_alpha(0.2)

    # add titling
    mainTitle = f"Mid-Atlantic Jan/Feb Snow Correlated to 500mb Heights (Lagged {lagtime * -1} Days)"
    subTitle = "\nData from ERA5 and xmACIS2 (Airports: IAD, BWI, PHL, EWR)"
    plt.title(mainTitle + subTitle, fontsize=9, weight='bold', loc='left')
    plt.title("DCAreaWx", fontsize=9, weight='bold', loc='right', color='gray')
    ax.text(-141.5, 1.5, 'Hatched regions indicate significant correlations', fontsize=9, weight='bold',
            transform=ccrs.PlateCarree())

    plt.savefig(rf"jan_snow_corr_{7 - lagtime * -1}.png", dpi=300, bbox_inches='tight')
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
    # plt.contourf(corr_data_array.longitude, corr_data_array.latitude, corr_data_array, 60, levels=np.arange(-0.15, 0.15, 0.005),
    #              extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
    # cbar = plt.colorbar(pad=0.017, aspect=25, shrink=0.84)
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
    # mainTitle = f"Mid-Atlantic Jan/Feb Snow Correlated to 500mb Heights (Lagged {lagtime * -1} Days)"
    # subTitle = "\nData from ERA5 and xmACIS2 (Airports: IAD, BWI, PHL, EWR)"
    # plt.title(mainTitle + subTitle, fontsize=9, weight='bold', loc='left')
    # plt.title("DCAreaWx", fontsize=9, weight='bold', loc='right', color='gray')
    # ax.text(182, 2, 'Hatched regions indicate significant correlations', fontsize=9, weight='bold',
    #         transform=ccrs.PlateCarree())
    #
    # plt.savefig(rf"jan_snow_corr_{7 - lagtime * -1}.png", dpi=300, bbox_inches='tight')
    # plt.show()


for i in range(-7, 1):
    plot(i, ['iad_snowfall_daily.csv', 'bwi_snowfall_daily.csv', 'phl_snowfall_daily.csv', 'ewr_snowfall_daily.csv'])
