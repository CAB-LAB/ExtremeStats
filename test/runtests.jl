using ExtremeStats
using Base.Test

# Test getSeasStat function
testx=[sin(linspace(0,2*pi,46)),2*sin(linspace(0,2*pi,46))]
msc,stdmsc=getSeasStat(testx,2,46)
@test all(msc.==1.5*sin(linspace(0,2*pi,46)))
@test all([isapprox(abs(stdmsc[i]),abs(0.5*sin(linspace(0,2*pi,46))[i])) for i = 1:46])
