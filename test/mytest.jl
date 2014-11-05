using ExtremeStats
println("Reading data")
@time x,lons,lats,nlon,nlat,NpY,N_years = load_X("/Net/Groups/BGI/people/mjung/FLUXCOM/_ENSEMBLES/8daily/TERall/","TERall_","TERall",2002:2007,[10,20],[55, 45]);
println("Calculating anomalies")
@time x=get_anomalies(x,NpY,nlon,nlat)
println("Labelling extremes")
@time el=label_Extremes(x,0.05);
println("Calculating time bounds")
@time getTbounds(el);
println("Calculating features")
@time f1=getFeatures(el,Features.Mean,Features.Max_z,Features.Min_z,Features.Duration,Features.Size,Features.NumPixel,Features.Quantile);
@time f2=getFeatures(el,Features.TS_ZValue,Features.TS_Area)
