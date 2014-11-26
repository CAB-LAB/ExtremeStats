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

function anomalies(series::Array,ilon,ilat,msc,stdmsc,smsc,sstdmsc,n,series_anomaly;NpY::Int=46,filt=9)
  # length of time series
  N    = size(series,3)
  @assert mod(N,NpY) == 0
  # divide length of time series by samples/year -> nr years
  NY=ifloor(N/NpY)
  # mean seasonal cycle
    getSeasStat(series,ilon,ilat,NY,NpY,msc,stdmsc,n)
  # smooth MAC
    smooth_circular!(msc, smsc, 9)
    #smooth_circular!(stdmsc, sstdmsc, 9)
  # deseasoanlized time series, difference of time series and mac
    for iyear=0:(NY-1),iday=1:NpY
        series_anomaly[ilon,ilat,iyear*NpY+iday]  = series[ilon,ilat,iyear*NpY+iday] - smsc[iday]
    end
end

function get_anomalies(x,NpY)
    anomar = similar(x)
    get_anomalies!(x,NpY,anomar)
    return(anomar)
end

function get_anomalies!(x,NpY)
    get_anomalies!(x,NpY,x)
    return(x)
end

function get_anomalies!(x,NpY,anomar)
    nlon=size(x,1)
    nlat=size(x,2)
    msc     = Array(Float64,NpY) #Allocate once
    stdmsc  = Array(Float64,NpY) #Allocate once
    smsc    = Array(Float64,NpY) #Allocate once
    sstdmsc = Array(Float64,NpY) #Allocate once
    n       = Array(Int64,NpY)
    for lon in 1:nlon, lat in 1:nlat
        anomalies(x,lon,lat,msc,stdmsc,smsc,sstdmsc,n,anomar,NpY=NpY)
        #println("$lon $lat")
    end
    return(anomar)
end
