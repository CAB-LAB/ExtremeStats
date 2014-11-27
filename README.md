# ExtremeStats

[![Build Status](https://travis-ci.org/meggart/ExtremeStats.jl.svg?branch=master)](https://travis-ci.org/meggart/ExtremeStats.jl)
This package is a collection of Extreme detection ans analysis algorithms implemented efficiently in the julia programming language.
Basic usage is as follows: 

Given a spatiotemporal data array **x** with dimensions in the order (lon,lat,time) you can call

    el=label_Extremes(x,kwargs...)

The following keyword arguments are valid:

* `quantile`: Quantile threshold to determine, which values are counted as extreme (tres=0.95 means the highest 5% values are considered extreme)
* `mask`: An additional mask having the same size as `x`, where only `true` values are considered for labelling
* `area`: A vector of the length `size(x,2)` giving the grid cell area depending on altitude. Defaults to equal area (one)
* `lons`: longitude valuse of the grid (optional)
* `lats`: lagtitude values of the grid (optional)
* `circular`: should th array be treated as closed along the longitude dimension as it is the case when analyzing the whole globe. defaults to true 

Returned is a structure of type `ExtremeList` containing a list with all the detected components.
After that you can use the getfeatures function to calculate some
statistics within each extreme

    getFeatures(el, features...)

where **el** is an ExtremeList followed by an optional number of features to be calculated. Currently implemented are the following features, please create a Gitlab issue if you wish more:

* Features.Mean
* Features.Max_z
* Features.Min_z
* Features.Duration
* Features.Size
* Features.NumPixel
* Features.Quantile 
* Features.TS_ZValue
* Features.TS_Area
* Features.Min_t
* Features.Max_t
* Features.Min_lon
* Features.Max_lon
* Features.Min_lat
* Features.Max_lat

This returns a Dataframe-like object with the requested features. 

In addition you have the possibility to subtract a smoothed mean annual cycle from the dataset using

    x2=get_anomalies(x,NpY)
    
This returns a new 3D array with only anomalies. There is also a mutating version of the function which does the anomaly-subtraction in-place.

    get_anomalies!(x,NpY)

Currently there is also a convenience function exported to load a BGI-like dataset. This is going to be replaced soon by using the MDIDS package. 

    load_X(data_path,fileprefix,varname,years,lon_range,lat_range)

This will read a datacube from NetCDF files in `data_path` with the given variable name for the given years, lons and lats. An example use would be:

    load_X("/Net/Groups/BGI/people/mjung/FLUXCOM/_ENSEMBLES/8daily/TERall/","TERall_","TERall",2001:2012,[20,40],[40, 20]);

If you want to combine a list of detected extremes into a big one, you can call

    elCombined = combineExtremes(el)

which returns a new ExtremeList having only a single Extreme. 
## Author(s)
Fabian Gans (BGI department)

## Credits
This package is heavily influenced by R code written by Sebastian Sippel. The anomaly detection algorithms are after a function written by Miguel Mahecha. 