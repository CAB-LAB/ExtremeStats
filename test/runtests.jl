using ExtremeStats
using Base.Test

# Test getSeasStat function
# Two subsequent years with 1st year sin(t) and 2nd year 2*sin(t)
t=linspace(0,2*pi,46)
testx=[sin(t),2*sin(t)]
testx=reshape(testx,1,1,92)
msc=Array(Float64,46)
stdmsc=Array(Float64,46)
n=Array(Int64,46)
ExtremeStats.getSeasStat(testx,1,1,2,46,msc,stdmsc,n)
# msc should be 1.5*si(t)
@test all(msc.==1.5*sin(t))
# std of msc should be 0.5 sin(t)
@test all([isapprox(abs(stdmsc[i]),abs(0.5*sin(t)[i])) for i = 1:46])
#Another test, direct comparison to R results
rtab=readdlm(joinpath(Pkg.dir(),"ExtremeStats/test/anomaly_example.csv"),' ')
r_in=reshape(rtab[:,1],1,1,552)
r_out=reshape(rtab[:,2],1,1,552)
j_out=get_anomalies(r_in,46)
@test all([isapprox(r_out[i],j_out[i]) for i=1:552])

#Test Extreme Cluster detection
a=zeros(Float32,100,100,100)
#Define an area vector
area=Float32[1:100]
#Define 2 clusters of very high values
b=Float32[i*0.1+k for i=1:10,j=1:10,k=1:10]
bweights=[area[j] for i=1:10,j=46:55,k=1:10]
a[91:100,46:55,66:75]=b
a[1:10,46:55,66:75]=b
#Also add a low extreme
c=rand(Float32,10,10,10)-100.0
a[51:60,51:60,51:60]=c
cweights=[area[j] for i=51:60,j=51:60,k=51:60]

#And add some NaNs
a[9:10,9:10,9:10]=NaN32
#Test without circular
el1=label_Extremes(a,quantile=0.95,area=area,circular=false)
f1=f1=getFeatures(el1,Features.Mean,Features.Max_z,Features.Min_z,Features.Duration,Features.Size,Features.NumPixel,
    Features.Min_t,Features.Max_t,Features.Min_lon,Features.Max_lon,Features.Min_lat,Features.Max_lat,Features.Quantile,
    Features.TS_ZValue, Features.TS_Area);
    
    #Now test all results
@test length(el1.extremes)==2
@test length(el1.extremes[1].zvalues)==length(b)
@test all(sort(el1.extremes[1].zvalues).==sort(b[:]))
@test all(sort(el1.extremes[2].zvalues).==sort(b[:]))
@test el1.extremes[1].tbounds==(66,75)
@test el1.extremes[1].lonbounds==(1,10)
@test el1.extremes[1].latbounds==(46,55)
@test el1.extremes[2].tbounds==(66,75)
@test el1.extremes[2].lonbounds==(91,100)
@test el1.extremes[2].latbounds==(46,55)
#Test first cluster
@test isapprox(f1[1][1],mean(b))
@test isapprox(f1[2][1],maximum(b))
@test isapprox(f1[3][1],minimum(b))
@test isapprox(f1[4][1],9)
@test isapprox(f1[5][1],mean(area[46:55])*1000)
@test f1[6][1]==1000
@test f1[7][1]==66
@test f1[8][1]==75
@test f1[9][1]==1
@test f1[10][1]==10
@test f1[11][1]==46
@test f1[12][1]==55
@test isapprox(f1[13][1][1],quantile(b[:],0.1),rtol=0.05)
@test isapprox(f1[13][1][2],quantile(b[:],0.5),rtol=0.05)
@test isapprox(f1[13][1][3],quantile(b[:],0.9),rtol=0.05)
@test all([isapprox(f1[14][1][i],(sum((b.*bweights),(1,2))./sum(bweights,(1,2)))[i]) for i=1:length(f1[14][1])])
@test all([isapprox(f1[15][1][i],sum(bweights,(1,2))[i]) for i=1:length(f1[15][1])])
#Test second cluster
@test isapprox(f1[1][2],mean(b))
@test isapprox(f1[2][2],maximum(b))
@test isapprox(f1[3][2],minimum(b))
@test isapprox(f1[4][2],9)
@test isapprox(f1[5][2],mean(area[46:55])*1000)
@test f1[6][2]==1000
@test f1[7][2]==66
@test f1[8][2]==75
@test f1[9][2]==91
@test f1[10][2]==100
@test f1[11][2]==46
@test f1[12][2]==55
@test isapprox(f1[13][2][1],quantile(b[:],0.1),rtol=0.05)
@test isapprox(f1[13][2][2],quantile(b[:],0.5),rtol=0.05)
@test isapprox(f1[13][2][3],quantile(b[:],0.9),rtol=0.05)
@test all([isapprox(f1[14][2][i],(sum((b.*bweights),(1,2))./sum(bweights,(1,2)))[i]) for i=1:length(f1[14][2])])
@test all([isapprox(f1[15][2][i],sum(bweights,(1,2))[i]) for i=1:length(f1[15][2])])



#Test with circular
el2 = label_Extremes(a,quantile=0.95,area=area,circular=true)
f2=getFeatures(el2,Features.Mean,Features.Max_z,Features.Min_z,Features.Duration,Features.Size,Features.NumPixel,
    Features.Min_t,Features.Max_t,Features.Min_lon,Features.Max_lon,Features.Min_lat,Features.Max_lat,Features.Quantile,
    Features.TS_ZValue, Features.TS_Area);
    
b=cat(1,b,b)
bweights=cat(1,bweights,bweights)
@test length(el2.extremes)==1
@test all(sort(el2.extremes[1].zvalues).==sort(b[:]))
@test isapprox(f2[1][1],mean(b))
@test isapprox(f2[2][1],maximum(b))
@test isapprox(f2[3][1],minimum(b))
@test isapprox(f2[4][1],9)
@test isapprox(f2[5][1],mean(area[46:55])*2000)
@test f2[6][1]==2000
@test f2[7][1]==66
@test f2[8][1]==75
@test f2[9][1]==1
@test f2[10][1]==100
@test f2[11][1]==46
@test f2[12][1]==55
@test isapprox(f2[13][1][1],quantile(b[:],0.1),rtol=0.05)
@test isapprox(f2[13][1][2],quantile(b[:],0.5),rtol=0.05)
@test isapprox(f2[13][1][3],quantile(b[:],0.9),rtol=0.05)
@test all([isapprox(f2[14][1][i],(sum((b.*bweights),(1,2))./sum(bweights,(1,2)))[i]) for i=1:length(f2[14][1])])
@test all([isapprox(f2[15][1][i],sum(bweights,(1,2))[i]) for i=1:length(f2[15][1])])

# Now label low extremes
el3=label_Extremes(a,quantile=0.05,area=area,circular=false)
f3=getFeatures(el3,Features.Mean,Features.Max_z,Features.Min_z,Features.Duration,Features.Size,Features.NumPixel,
    Features.Min_t,Features.Max_t,Features.Min_lon,Features.Max_lon,Features.Min_lat,Features.Max_lat,Features.Quantile,
    Features.TS_ZValue, Features.TS_Area);
@test length(el3.extremes)==1
@test all(sort(el3.extremes[1].zvalues).==sort(c[:]))
@test isapprox(f3[1][1],mean(c))
@test isapprox(f3[2][1],maximum(c))
@test isapprox(f3[3][1],minimum(c))
@test isapprox(f3[4][1],9)
@test isapprox(f3[5][1],mean(area[51:60])*1000)
@test f3[6][1]==1000
@test f3[7][1]==51
@test f3[8][1]==60
@test f3[9][1]==51
@test f3[10][1]==60
@test f3[11][1]==51
@test f3[12][1]==60
@test isapprox(f3[13][1][1],quantile(c[:],0.1),rtol=0.05)
@test isapprox(f3[13][1][2],quantile(c[:],0.5),rtol=0.05)
@test isapprox(f3[13][1][3],quantile(c[:],0.9),rtol=0.05)
@test all([isapprox(f3[14][1][i],(sum((c.*cweights),(1,2))./sum(cweights,(1,2)))[i]) for i=1:length(f3[14][1])])
@test all([isapprox(f3[15][1][i],sum(cweights,(1,2))[i]) for i=1:length(f3[15][1])])
