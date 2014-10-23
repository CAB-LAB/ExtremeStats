module ExtremeStats
export anomalies, get_anomalies, Extreme, load_X, label_Extremes, ExtremeList
import Images.label_components
import NetCDF.ncread
include("Features.jl")

function load_X(evapofrac_path,years,lat_range,lon_range)

  N_years   = length(years)
  y=years[1]
  lons=ncread(joinpath(evapofrac_path,"EvapoFrac.$(years[1]).nc"),"lon")
  lats=ncread(joinpath(evapofrac_path,"EvapoFrac.$(years[1]).nc"),"lat")
  time=ncread(joinpath(evapofrac_path,"EvapoFrac.$(years[1]).nc"),"time")
  NpY=length(time)
  lonstep=lons[2]-lons[1]
  latstep=lats[2]-lats[1]
  cmplon=lonstep>0 ? .> : .<
  cmplat=latstep>0 ? .> : .<
  ilon_range=sort!([findfirst(cmplon(lons,lon_range[1])),findfirst(cmplon(lons,lon_range[2]-lonstep))])
  ilat_range=sort!([findfirst(cmplat(lats,lat_range[1])),findfirst(cmplat(lats,lat_range[2]-latstep))])
  nlon=ilon_range[2]-ilon_range[1]
  nlat=ilat_range[2]-ilat_range[1]
  lons=lons[ilon_range[1]:ilon_range[2]];
  lats=lats[ilat_range[1]:ilat_range[2]];
  x=Array(Float32,nlon,nlat,N_years*46)
  iyear=1
  for iyear in 1:N_years
      println("$iyear")
      x[:,:,((iyear-1)*NpY+1):iyear*NpY]=ncread(joinpath(evapofrac_path,"EvapoFrac.$(years[iyear]).nc"),"EvapoFrac",[ilon_range[1],ilat_range[1],1],[nlon,nlat,-1]);
  end
  return(x,lons,lats,nlon,nlat,NpY,N_years)
end

random_x(nlon=200,nlat=200,N_years=7)=rand(Float32,nlon,nlat,N_years*46);


function getSeasStat(series,NY,NpY,msc,stdmsc)
    NY*NpY==length(series) || error("Length of series not good")
    for iday = 1:NpY
        s=zero(eltype(series))
        s2=zero(eltype(series))
        n=0
        for iyear=0:(NY-1)
            v=series[(iyear)*NpY+iday]
            if !isnan(v)
                s2=s2+v*v
                s=s+v
                n=n+1
            end
        end
        msc[iday]=s/n
        stdmsc[iday]=sqrt(s2/n-s*s/(n*n))
    end
    msc,stdmsc
end

function smooth_circular!(xin::Vector,xout::Vector,wl::Integer=9)
    mod(wl,2)==1 || error("Window length must be uneven")
    offs = div(wl,2)
    fwl=float64(wl)
    n=length(xin)
    awl=0
    s=0.0
    for i=(n-offs+1):n
        if !isnan(xin[i])
            s=s+xin[i]
            awl=awl+1
        end
    end
    for i=1:(offs+1)
        if !isnan(xin[i])
            s=s+xin[i];
            awl=awl+1
        end
    end
    #Start main loop
    for i=1:offs
        xout[i]=s/awl
        if !isnan(xin[n-offs+i])
            s=s-xin[n-offs+i]
            awl=awl-1
        end
        if !isnan(xin[offs+1+i])
            s=s+xin[offs+1+i]
            awl=awl+1
        end
    end
    for i=(offs+1):(n-offs-1)
        xout[i]=s/awl
        if !isnan(xin[i-offs])
            s=s-xin[i-offs]
            awl=awl-1
        end
        if !isnan(xin[offs+1+i])
            s=s+xin[offs+1+i]
            awl=awl+1
        end
    end
    for i=(n-offs):n
        xout[i]=s/awl
        if !isnan(xin[i-offs])
            s=s-xin[i-offs]
            awl=awl-1
        end
        if !isnan(xin[i-n+offs+1])
            s=s+xin[i-n+offs+1]
            awl=awl+1
        end
    end
end

function anomalies(series::SubArray,msc,stdmsc,smsc,sstdmsc,series_anomaly;NpY::Int=46,filt=9)
  # length of time series
  N    = length(series)
  @assert mod(N,NpY)==0
  # devide length of time series by samples/year -> nr years
  NY=ifloor(N/NpY)
  # mean seasonal cycle
    getSeasStat(series,NY,NpY,msc,stdmsc)
  # smooth MAC
    smooth_circular!(msc, smsc, 9)
    smooth_circular!(stdmsc, sstdmsc, 9)
  # deseasoanlized time series, difference of time series and mac
    for iyear=1:NY,iday=1:NpY
        series_anomaly[(iyear-1)*NpY+iday]  = (series[(iyear-1)*NpY+iday] - smsc[iday])/sstdmsc[iday]
    end
end

function nanquantile(xtest,q)
nNaN  = mapreduce(isnan,+,xtest)
i=iround((length(xtest)-nNaN)*q)
select(xtest[:],i)
end

function label_Extremes{T}(x::Array{T,3},quantile::Number;circular::Bool=false,pattern::BitArray=trues(3,3,3))
  nlon=size(x,1)
  #First calculate threshold
  cmpfun = (quantile < 0.5) ? .< : .>
  tres = nanquantile(x,quantile)
  #Then allocate bitArray that hold the true/falses
  offs = circular ? 1 : 0
  x2=BitArray(size(x,1)+offs,size(x,2),size(x,3))
  #Fill BitArray
  x2[1:nlon,:,:]=cmpfun(x,tres)
  if circular
    #Now attach first slice to the end (circular globe)
    x2[nlon+1,:,:]=x2[1,:,:];
  end
  #Do the labelling
  lx=label_components(x2,pattern);
  x2=0;gc();
  circular && renameLabels(lx,x);
  el=ExtremeList(x,lx);
  return(el)
end



function get_anomalies(x,NpY,nlon,nlat)
msc    = Array(Float64,NpY) #Allocate once
stdmsc = Array(Float64,NpY) #Allocate once
smsc    = Array(Float64,NpY) #Allocate once
sstdmsc = Array(Float64,NpY) #Allocate once
anomar=similar(x)
for lon in 1:nlon, lat in 1:nlat
    anomalies(sub(x,lon,lat,:),msc,stdmsc,smsc,sstdmsc,sub(anomar,lon,lat,:),NpY=NpY)
    #println("$lon $lat")
end
  return(anomar)
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

type Extreme{T}
  index::Int64
  locs::Array{Int,2}
  zvalues::Array{T,1}
end



type ExtremeList{T}
  extremes::Array{Extreme{T}}
  features::Dict{Symbol, Array{Float32,1}}
end

typealias FeatureVector{T} Vector{T}
function ExtremeList(x,lx)
  nEx=maximum(lx)
  numCells=countNumCell(lx,nEx)
  nempty=sum(numCells.==0);
  extremeList=[Extreme(i,zeros(Int,numCells[i],3),Array(eltype(x),numCells[i])) for i=1:nEx];
  indices2List(lx,x,extremeList);
  o=sortperm(numCells,rev=true)
  extremeList=extremeList[o];
  deleteat!(extremeList,(nEx-nempty+1):nEx);
  return(ExtremeList(extremeList,Dict{Symbol, FeatureVector{eltype(x)}}()))
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


end # module

