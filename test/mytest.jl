using ExtremeStats

x,lons,lats,nlon,nlat,NpY,N_years = load_X("/Net/Groups/BGI/people/uweber/_data/2MM/EvapoFrac_FLuxcom/Data/",2002:2007,[45, 55],[0, 20]);
el=label_Extremes(x,0.05);
ExtremeStats.defineFeatures(Features.Mean,Features.Max_z,Features.Min_z)
@time getFeatures(el)
