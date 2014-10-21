using ExtremeStats
using Images
NpY=46
nlon=200
nlat=200
NY=7
x=ExtremeStats.random_x(nlon,nlat,NY)
inan=rand(1:(nlon*nlat*NY),div(nlon*nlat*NY,20))
x[inan]=NaN


@time anomar=get_anomalies(x,NpY,nlon,nlat)
@time tres=nanquantile(x,0.95)
x2=x.>tres
@time lx=label_components(x2,trues(3,3,3))
minimum(lx)
nEx=maximum(lx)
nFeat=4
numCells=countNumCell(lx,nEx)
@time extremeList=[Extreme(i,Array(Int,3,numCells[i]),Array(Float64,nFeat)) for i=1:nEx]

function indices2List(larr,extremeList)
  curind=ones(length(larr))
  z=1;
  for (i,j,k) in (1:size(larr,1), 1:size(larr,2), 1:size(larr,3))
    if larr[i,j,k]>0
      extremeList[z].locs[1,curind[larr[i,j,k]]]=i
      extremeList[z].locs[2,curind[larr[i,j,k]]]=j
      extremeList[z].locs[3,curind[larr[i,j,k]]]=k
    end
    z=z+1
    curind[]
  end
end
indices2List(lx,extremeList)
extremeList[500]
extremeList=extremeList[sortperm(numCells)]
extremeList[end-1]
