# ExtremeStats

[![Build Status](https://travis-ci.org/meggart/ExtremeStats.jl.svg?branch=master)](https://travis-ci.org/meggart/ExtremeStats.jl)
This package is a collection of Extreme detection ans analysis algorithms implemented efficiently in the julia programming language.
Basic usage is as follows: 

Given a spatiotemporal data array **x** with dimensions in the order (lon,lat,time) you can call

    el=label_Extremes(x,tres,area=area,lons=lons,lats=lats,circular=true)

where **tres** is a threshold quantile to determine, which values are extreme (tres=0.95 means the highest 5% values are considered extreme). 
You can also provide the grid cell area as a vector of length nlatitudes. The **circular** argument determines, if the array should be treated as 
closed along the longitude dimension as it is the case when analyzing the whole globe. 

Returned is a structure of type ExtremeList, containing a Vector of all connected Extremes sorted by their number of pixels. You should call the function

    getTbounds(el)

to calculate 3D bounding boxes for each extreme prior to calculating features. After that you can use the getfeatures function to calculate some
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
    x2=get_anomalies(x,NpY,nlon,nlat)
This returns a new 3D array with only anomalies. 

## Author(s)
Fabian Gans (BGI department)

## Credits
This package is heavily influenced by R code written by Sebastian Sippel. The anomaly detection algorithms are after function written by Miguel Mahecha. 