using ExtremeStats
using Base.Test

# Test getSeasStat function
# Two subsequent years with 1st year sin(t) and 2nd year 2*sin(t)
t=linspace(0,2*pi,46)
testx=[sin(t),2*sin(t)]
msc=Array(Float64,46)
stdmsc=Array(Float64,46)
getSeasStat(testx,2,46,msc,stdmsc)
# msc should be 1.5*si(t)
@test all(msc.==1.5*sin(t))
# std of msc should be 0.5 sin(t)
@test all([isapprox(abs(stdmsc[i]),abs(0.5*sin(t)[i])) for i = 1:46])


#Test Extreme Cluster detection
a=zeros(Float32,10,10,10)
#Define 2 clusters of very high values
a[1:2,5:6,5:6]=100.0
a[9:10,5:6,5:6]=100.0
#And add some NaNs
a[9:10,9:10,9:10]=NaN32
#Test without circular
el1=label_Extremes(a,0.95,circular=false)
@test length(el1.extremes)==2
@test length(el1.extremes[1].zvalues)==8
@test all(el1.extremes[1].zvalues.==100.0)
@test length(el1.extremes[2].zvalues)==8
@test all(el1.extremes[2].zvalues.==100.0)

#Test with circular
el2=label_Extremes(a,0.95,circular=true)
@test length(el1.extremes)==1
@test length(el1.extremes[1].zvalues)==16
@test all(el1.extremes[1].zvalues.==100.0)

