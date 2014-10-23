using ExtremeStats
using Images
using PyPlot

NpY=10
nlon=20
nlat=20
NY=1
x=ExtremeStats.random_x(nlon,nlat,NY)
inan=rand(1:(nlon*nlat*NY),div(nlon*nlat*NY,20))
x[inan]=NaN
area=ones(Float32,nlat)



@time anomar=get_anomalies(x,NpY,nlon,nlat)
@time tres=nanquantile(x,0.95)
x2=BitArray(nlon+1,nlat,NY*NpY)
x2[1:nlon,:,:]=x.>tres
#Now attach first slice to the end (circular globe)
x2[nlon+1,:,:]=x2[1,:,:]

@time lx=label_components(x2,trues(3,3,3))
#Take the last slice and compare it to the first
renames=Dict{Int,Int}()
for k=1:nlat, t=1:(NY*NpY)
  if (lx[nlon+1,k,t]>0) && (!haskey(renames,lx[nlon+1,k,t]))
    renames[lx[nlon+1,k,t]]=lx[1,k,t]
  end
end
i=nlon+1;sren=1;
while sren>0
  sren=0;
  for k=1:nlat, l=1:(NY*NpY)
    if haskey(renames,lx[i,k,l])
      lx[i,k,l]=renames[lx[i,k,l]]
      sren=sren+1;
    end
  end
  i=i-1;
  println("$i $sren")
end
length(renames)
nEx=maximum(lx)
numCells=countNumCell(lx,nEx)
sum(numCells.==0)
inds=findnz(numCells)
@time extremeList=[Extreme(i,Array(Int,3,numCells[i]),Array(eltype(x),numCells[i])) for i=1:nEx]

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
