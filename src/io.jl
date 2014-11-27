import NetCDF.nccreate, NetCDF.ncwrite, NetCDF.ncclose
# writes the first 1000 extremes of an extremelist to a nectdf file
function writeExtremes(el::ExtremeList,folder::String;nmax=1000,prefix="e")
    run(`mkdir -p $(folder)`)
    for iex=1:min(nmax,length(el.extremes))
        si=@sprintf("%04d",iex)
        filename="$(folder)/$(prefix)$(si).nc"
        isfile(filename) && rm(filename)
        nccreate(filename,"locs",{"units"=>"Lon,Lat,Time"},"pixel",[1:size(el.extremes[iex].locs,1)],"dimension",[1,2,3,4])
        ncwrite(el.extremes[iex].locs,filename,"locs",start=[1,1],count=[-1,3])
        ncwrite(el.extremes[iex].zvalues,filename,"locs",start=[1,4],count=[-1,1])
        ncclose()
    end
end

function unroll_tuples(f::Vector)
    nwrites=Int[length(fi[1]) for fi in f]
    maximum(nwrites)>5 && error("Do not use the writeFeatures function to write time series, use it only for single indices")
    nEx=length(f[1])
    fnew=Array(Any,sum(nwrites))
    inew=1
    for i=1:length(f)
        l=length(f[i][1])
        if l==1
            #Everything is ok
            fnew[inew]=f[i]
            inew=inew+1
        else
            for j=1:l
                fnew[inew]=[f[i][k][j] for k=1:nEx]
                inew=inew+1
            end
        end
    end
    return fnew
end
                
    

function writeFeatures(f::Vector,featnames::Vector,file;overwrite=true)
    fnew=unroll_tuples(f)
    length(fnew)==length(featnames) || error("Number of feature names must match number of written features. Make sure that you named each quantile separately!")
    featDict=Dict{Any,Any}()
    [featDict["$i"]="$(featnames[i])" for i=1:length(featnames)]
    #Write features
    isfile(file) && overwrite && rm(file)
    nccreate(file,"features",{"unit"=>"diverse"},"N_Extreme",[1:length(fnew[1])],"Feat",[1:length(featnames)],featDict)
    for i=1:length(fnew)
            ncwrite(fnew[i],file,"features",start=[1,i],count=[-1,1])
    end
    ncclose()
end

writeTimeSeries(f::Vector,filename,varNames::String)=writeTimeSeries(f,filename,[varNames])
function writeTimeSeries(f::Vector,filename,varNames::Vector)
    isfile(filename) && rm(filename)
    ntstep=length(f[1][1])
    length(varNames)==length(f) || error("number of variable names must equal")
    for i=1:length(f)
        nccreate(filename,varNames[i],{"unit"=>"none"},"timestep",[1:ntstep])
    end
    for i=1:length(f)
        ncwrite(f[i][1],filename,varNames[i])
    end
    ncclose()
end