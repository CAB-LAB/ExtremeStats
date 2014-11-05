using ExtremeStats
println("Reading data")
@time x,lons,lats,nlon,nlat,NpY,N_years = load_X("/Net/Groups/BGI/people/mjung/FLUXCOM/_ENSEMBLES/8daily/TERall/","TERall_","TERall",2001:2012,[20,40],[40, 20]);
println("Calculating anomalies")
gc()
@time x=get_anomalies(x,NpY,nlon,nlat)
println("Labelling extremes")
gc()
@time el=label_Extremes(x,0.05);
println("Calculating time bounds")
gc()
@time getTbounds(el);
println("Calculating features")
@time f1=getFeatures(el,Features.Mean,Features.Max_z,Features.Min_z,Features.Duration,Features.Size,Features.NumPixel,Features.Quantile);
gc()
@time f2=getFeatures(el,Features.TS_ZValue,Features.TS_Area)
