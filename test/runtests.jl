using ExtremeStats
using Base.Test

# Test getSeasStat function
# Two subsequent years with 1st year sin(t) and 2nd year 2*sin(t)
t=linspace(0,2*pi,46)
testx=[sin(t),2*sin(t)]
testx=reshape(testx,1,1,92)
msc=Array(Float64,46)
stdmsc=Array(Float64,46)
getSeasStat(testx,1,1,2,46,msc,stdmsc)
# msc should be 1.5*si(t)
@test all(msc.==1.5*sin(t))
# std of msc should be 0.5 sin(t)
@test all([isapprox(abs(stdmsc[i]),abs(0.5*sin(t)[i])) for i = 1:46])


#Test Extreme Cluster detection
a=zeros(Float32,100,100,100)
#Define an area vector
area=Float32[1:100]
#Define 2 clusters of very high values
b=Float32[i*0.1+k for i=1:10,j=1:10,k=1:10]
a[91:100,46:55,66:75]=b
a[1:10,46:55,66:75]=b
#And add some NaNs
a[9:10,9:10,9:10]=NaN32
#Test without circular
el1=label_Extremes(a,0.95,area=area,circular=false)
getTbounds(el1);
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
@test isapprox(f1[1][1],mean(b))
@test isapprox(f1[2][1],maximum(b))
@test isapprox(f1[3][1],minimum(b))
@test isapprox(f1[4][1],9)


#Test with circular
el2=label_Extremes(a,0.95,circular=true)
@test length(el2.extremes)==1
@test all(sort(el2.extremes[1].zvalues).==sort([b[:], b[:]]))

