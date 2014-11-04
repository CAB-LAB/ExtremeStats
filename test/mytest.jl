using ExtremeStats

x,lons,lats,nlon,nlat,NpY,N_years = load_X("/Users/fgans/scratch/data/",2002:2007,[45,55],[25, 45]);
el=label_Extremes(x,0.05);
@time f=getFeatures(el,Features.Mean,Features.Max_z,Features.Min_z,Features.Duration,Features.Size,Features.NumPixel,Features.Quantile)
