module ExtremeStats

include("anomalies.jl")
include("features.jl")

export get_anomalies, get_anomalies!, Extreme, load_X, label_Extremes, ExtremeList, Features, getFeatures, combineExtremes, sortby, writeExtremes, writeFeatures, writeTimeSeries
import Images.label_components
import NetCDF.ncread

type Extreme{T}
  index::Int64
  locs::Array{Int,2}
  zvalues::Array{T,1}
  tbounds::(Int,Int)
  lonbounds::(Int,Int)
  latbounds::(Int,Int)
end



type ExtremeList{T,U,V}
  extremes::Array{Extreme{T}}
  area::Vector{U}
  lons::Vector{V}
  lats::Vector{V}
end

include("io.jl")

getindex(el::ExtremeList,i...)=getindex(el.extremes,i...)
function sortby!(el::ExtremeList,x::Vector;rev=false)
    o=sortperm(x,rev=rev)
    el.extremes[:]=el.extremes[o]
end


random_x(nlon=200,nlat=200,N_years=7)=rand(Float32,nlon,nlat,N_years*46);


function getSeasStat(series,ilon,ilat,NY,NpY,msc,stdmsc,n)
    fill!(msc,0.0)
    fill!(stdmsc,0.0)
    fill!(n,0.0)
    for iyear=0:(NY-1)
        for iday = 1:NpY
            v=series[ilon,ilat,iyear*NpY+iday]
            if !isnan(v)
                stdmsc[iday]=stdmsc[iday]+v*v
                msc[iday]=msc[iday]+v
                n[iday]=n[iday]+1
            end
        end
    end
    for iday = 1:NpY 
        if n[iday]>0
            stdmsc[iday]=sqrt(max(stdmsc[iday]/n[iday]-msc[iday]*msc[iday]/(n[iday]*n[iday]),zero(eltype(series))))
            msc[iday]=msc[iday]/n[iday]
        else
            msc[iday]=nan(eltype(series))
            stdmsc[iday]=nan(eltype(series))
        end
    end
    msc,stdmsc
end


function nanquantile(xtest,q;useabs=false) 
    x = Array(eltype(xtest),length(xtest))
    nanquantile!(xtest,q,x,useabs=useabs)
end

function nanquantile!(xtest,q,x;useabs=false) 
    useabs ? for i=1:length(xtest) x[i]=abs(xtest[i]) end : copy!(x,xtest)
    nNaN  = 0; for i=1:length(xtest) nNaN = isnan(xtest[i]) ? nNaN+1 : nNaN end
    lv=length(xtest)-nNaN
    lv==0 && return(nan(eltype(x)))
    index = 1 + (lv-1)*q
    lo = ifloor(index)
    hi = iceil(index)
    vals=select!(x,lo:hi)
    h=index - lo
    length(vals)==1 && return(vals[1])
    r = (1.0-h)*vals[1] + h*vals[2]
end


function label_Extremes{T}(x::Array{T,3};quantile::Number=0.0,mask::BitArray{3}=trues(1,1,1),circular::Bool=false,pattern::BitArray=trues(3,3,3),area=ones(Float32,size(x,2)),lons=linspace(0,360,size(x,1)),lats=linspace(90,-90,size(x,2)),normscope="global")
    
    size(x)==size(mask) || size(mask)==(1,1,1) || error("Sizes of x and mask must be the same")
    nlon,nlat,ntime = size(x)
    offs = circular ? 1 : 0
    
    if quantile==one(quantile) || quantile==zero(quantile)
        x2=trues(size(x,1)+offs,size(x,2),size(x,3))
    else
        cmpfun = (quantile < 0.5) ? .< : .>
        #Then allocate bitArray that hold the true/falses
        x2=BitArray(size(x,1)+offs,size(x,2),size(x,3))
        if normscope=="global"
          #First calculate threshold
          tres = nanquantile(x,quantile)
          #Fill BitArray
          x2[1:nlon,:,:]=cmpfun(x,tres)
        elseif normscope=="local"
          qar1      = zeros(eltype(x),ntime)
          qar2      = zeros(eltype(x),ntime)
          for ilat=1:nlat, ilon=1:nlon
            for itime = 1:ntime
              qar1[itime] = x[ilon,ilat,itime]
            end
            tres=nanquantile!(qar1,quantile,qar2)
            x2[ilon,ilat,:]=cmpfun(qar1,tres)
            
          end
        end
    end
    
    if size(mask) != (1,1,1)
        presum=sum(x2)
        for i=1:size(mask,3), j=1:size(mask,2), k=1:size(mask,1)
            x2[k,j,i]=x2[k,j,i] && mask[k,j,i]
        end
        aftersum=sum(x2)
        println("$(presum-aftersum) of $(presum) elements were removed by additional condition ()$((presum-aftersum)/presum)%)")
    end

    if circular
      #Now attach first slice to the end (circular globe)
      x2[nlon+1,:,:]=x2[1,:,:];
    end
  #Do the labelling
  lx=label_components(x2,pattern);
  x2=0;gc();
  circular && renameLabels(lx,x);
  el=ExtremeList(x,lx,area,lons,lats);
  return(el)
end


function indices2List(larr,xarr,extremeList)
    curind=ones(Int,length(extremeList))
    for i=1:size(xarr,3), j=1:size(xarr,2), k=1:size(xarr,1)
        if larr[k,j,i]>0
            extremeList[larr[k,j,i]].locs[curind[larr[k,j,i]],1]=k
            extremeList[larr[k,j,i]].locs[curind[larr[k,j,i]],2]=j
            extremeList[larr[k,j,i]].locs[curind[larr[k,j,i]],3]=i
            extremeList[larr[k,j,i]].zvalues[curind[larr[k,j,i]]]=xarr[k,j,i]
            curind[larr[k,j,i]]=curind[larr[k,j,i]]+1
    end
  end
    return(curind)
end


function getTbounds(el::ExtremeList)
  for e in el.extremes
    tmin=typemax(Int)
    tmax=0
    lonmin=typemax(Int)
    lonmax=0
    latmin=typemax(Int)
    latmax=0
    for i=1:length(e.zvalues)
      if e.locs[i,3]<tmin tmin=e.locs[i,3] end
      if e.locs[i,3]>tmax tmax=e.locs[i,3] end
      if e.locs[i,1]<lonmin lonmin=e.locs[i,1] end
      if e.locs[i,1]>lonmax lonmax=e.locs[i,1] end
      if e.locs[i,2]<latmin latmin=e.locs[i,2] end
      if e.locs[i,2]>latmax latmax=e.locs[i,2] end
    end
    e.tbounds=(tmin,tmax)
    e.lonbounds=(lonmin,lonmax)
    e.latbounds=(latmin,latmax)
  end
end

typealias FeatureVector{T} Vector{T}

function ExtremeList{T}(x::Array{T,3},lx::Array{Int,3},area=ones(Float32,size(x,2)),lons=linspace(0,360,size(x,1)),lats=linspace(90,-90,size(x,2)))
  nEx=maximum(lx)
  numCells=countNumCell(lx,nEx)
  nempty=sum(numCells.==0);
  extremeList=[Extreme(i,zeros(Int,numCells[i],3),Array(eltype(x),numCells[i]),(0,0),(0,0),(0,0)) for i=1:nEx];
  indices2List(lx,x,extremeList);
  o=sortperm(numCells,rev=true)
  extremeList=extremeList[o];
  deleteat!(extremeList,(nEx-nempty+1):nEx);
  return(ExtremeList(extremeList,area,lons,lats))
end

function combineExtremes(elin::ExtremeList;nEx=length(elin.extremes))
    totlen=0
    for i=1:nEx
        totlen=totlen+length(elin.extremes[i].zvalues)
    end
    locs=Array(Int,totlen,3)
    zvalues=Array(eltype(elin.extremes[1].zvalues),totlen)
    k=1
    for i=1:nEx
        for j=1:length(elin.extremes[i].zvalues)
            locs[k,1]=elin.extremes[i].locs[j,1]
            locs[k,2]=elin.extremes[i].locs[j,2]
            locs[k,3]=elin.extremes[i].locs[j,3]
            zvalues[k]=elin.extremes[i].zvalues[j]
            k=k+1
        end
    end
    mins=minimum(locs,1)
    maxs=maximum(locs,1)
    e=Extreme(1,locs,zvalues,(mins[3],maxs[3]),(mins[1],maxs[1]),(mins[2],maxs[2]))
    elout=ExtremeList([e],elin.area,elin.lons,elin.lats)
end

function countNumCell(labelList,nEx)
  lAr=zeros(Int,nEx)
  for i=1:length(labelList)
    j=labelList[i]
    if j>0
      lAr[j]=lAr[j]+1
    end
  end
  return lAr
end

function renameLabels(lx,x)
  #Check if we have to relabel, if longitudes are padded
  nlon=size(x,1)
  nlat=size(x,2)
  ntime=size(x,3)
  size(lx,1)==size(x,1) && return(lx)
  size(lx,1)==(nlon+1) || error("Something is wrong with the lon dimensions")
  size(lx,2)==nlat || error("Something is wrong with the lat dimensions")
  size(lx,3)==ntime || error("Something is wrong with the time dimensions")
  renames=Dict{Int,Int}()
  for k=1:nlat, t=1:ntime
    if (lx[nlon+1,k,t]>0) && (!haskey(renames,lx[nlon+1,k,t]))
      renames[lx[nlon+1,k,t]]=lx[1,k,t]
    end
  end
  i=nlon+1;sren=1;
  while sren>0 && i>0
    sren=0;
    for k=1:nlat, l=1:(ntime)
      if haskey(renames,lx[i,k,l])
        lx[i,k,l]=renames[lx[i,k,l]]
        sren=sren+1;
      end
    end
    i=i-1;
  end
  lx[nlon+1,:,:]=0
end


function getFeatures(el::ExtremeList,features...)
    el.extremes[1].tbounds==(0,0) && getTbounds(el)
    myf       = ExtremeStats.Features.calcFeatureFunction(features...)
    prearrays = ExtremeStats.Features.getPreArrays(length(el.extremes[1].zvalues),eltype(el.extremes[1].zvalues),features...)
    eval(myf)
    retar     = [Array(Features.rettype(f,el),length(el.extremes)) for f in features]
    for i=1:length(el.extremes)
        ret=getFeatures(el.extremes[i],el.area,el.lons,el.lats,prearrays)
        for j=1:length(ret) retar[j][i]=ret[j] end
    end
    return retar
end

end # module

