module Features
const chosenFeatures=Array(DataType,0)
#Define Accumulator behaviors
abstract Accumulator
init(t::Type{Accumulator}) = error("Init behavior for Accumulator $t not defined.")
accu(t::Type{Accumulator}) = error("Accumulation behavior for Accumulator $t not defined.")
type Max_z_acc <: Accumulator end
init(::Type{Max_z_acc})= :(max_z_acc=typemin(eltype(z)))
accu(::Type{Max_z_acc})= :(if z[i]>max_z_acc max_z_acc=z[i] end)
type Min_z_acc <: Accumulator end
init(::Type{Min_z_acc})= :(min_z_acc=typemax(eltype(z)))
accu(::Type{Min_z_acc})= :(if z[i]<min_z_acc min_z_acc=z[i] end)
type Sum_z_acc <: Accumulator end
init(::Type{Sum_z_acc})= :(sum_z_acc=zero(eltype(z)))
accu(::Type{Sum_z_acc})= :(sum_z_acc=sum_z_acc+area[x[i,2]]*z[i])
type Size_acc <: Accumulator end
init(::Type{Size_acc})= :(size_acc=zero(eltype(area)))
accu(::Type{Size_acc})= :(size_acc=size_acc+area[x[i,2]])
type Weight_acc <: Accumulator end
init(::Type{Weight_acc})= :(weight_acc=zero(eltype(area)))
accu(::Type{Weight_acc})= :(weight_acc=weight_acc+area[x[i,2]])
#Time series accumulators are not preallocated and so the performance bottleneck at the moment....
type Z_TS_acc <: Accumulator end
init(::Type{Z_TS_acc})= :(z_ts_acc=zeros(eltype(z),max_t-min_t+1))
accu(::Type{Z_TS_acc})= :(z_ts_acc[x[i,3]-min_t+1]=z_ts_acc[x[i,3]-min_t+1]+z[i]*area[x[i,2]])
type AR_TS_acc <: Accumulator end
init(::Type{AR_TS_acc})= :(ar_ts_acc=zeros(eltype(area),max_t-min_t+1))
accu(::Type{AR_TS_acc})= :(ar_ts_acc[x[i,3]-min_t+1]=ar_ts_acc[x[i,3]-min_t+1]+area[x[i,2]])


type Mean end
accumulators(::Type{Mean})=(Sum_z_acc,Weight_acc)
final(::Type{Mean})= :(sum_z_acc/weight_acc)
rettype(::Type{Mean},e)=eltype(e.extremes[1].zvalues)
type Max_z end
accumulators(::Type{Max_z})=(Max_z_acc,)
final(::Type{Max_z})= :(max_z_acc)
rettype(::Type{Max_z},e)=eltype(e.extremes[1].zvalues)
type Min_z end
accumulators(::Type{Min_z})=(Min_z_acc,)
final(::Type{Min_z})= :(min_z_acc)
rettype(::Type{Min_z},e)=eltype(e.extremes[1].zvalues)
type Duration end
accumulators(::Type{Duration})=()
final(::Type{Duration})= :(max_t-min_t)
rettype(::Type{Duration},e)=Int
type Size end
accumulators(::Type{Size})=(Size_acc,)
final(::Type{Size})= :(size_acc)
rettype(::Type{Size},e)=eltype(e.area)
type NumPixel end
accumulators(::Type{NumPixel})=()
final(::Type{NumPixel})= :(length(z))
rettype(::Type{NumPixel},e)=Int
type TS_Area end
accumulators(::Type{TS_Area})=(AR_TS_acc,)
final(::Type{TS_Area})= :(ar_ts_acc)
rettype(::Type{TS_Area},e)=Vector{eltype(e.area)}
type TS_ZValue end
accumulators(::Type{TS_ZValue})=(AR_TS_acc,Z_TS_acc)
final(::Type{TS_ZValue})= quote
  begin
    for i_ts_z=1:length(ar_ts_acc)
      z_ts_acc[i_ts_z]=z_ts_acc[i_ts_z]/ar_ts_acc[i_ts_z]
    end
    z_ts_acc
  end
end
rettype(::Type{TS_ZValue},e)=Vector{eltype(e.extremes[1].zvalues)}


#Now the fun begins! quantile calculations!
#We use preallocated arrays to be fast and choose one of the according size
type Quantile end
accumulators(::Type{Quantile})=(Max_z_acc,)
final(::Type{Quantile})=quote
  begin
    l=length(z)
    if l==1
      (z[1],z[1],z[1])
    else
      imagnitude=iceil(log10(l))
      ar=prearrays[imagnitude]
      fill!(ar,max_z_acc)
      copy!(ar,z)
      sort!(ar)
      i1=min(iround(l*0.1),l)
      i1=max(i1,1)
      i2=min(iround(l*0.5),l)
      i2=max(i2,1)
      i3=min(iround(l*0.9),l)
      i3=max(i3,1)
      (ar[i1],ar[i2],ar[i3])
    end
  end
end
rettype(::Type{Quantile},e)=(eltype(e.extremes[1].zvalues),eltype(e.extremes[1].zvalues),eltype(e.extremes[1].zvalues))

function calcFeatureFunction(features::DataType...)
    # Get required accumulators
    empty!(chosenFeatures)
    accus = Array(Any,0)
    for feat in features
        push!(chosenFeatures,feat)
        for a in accumulators(feat)
            in(a,accus) || push!(accus,a)
        end
    end
    functionbody = quote
        x=e.locs
        z=e.zvalues
        n=length(z)
        (min_t,max_t)=e.tbounds
    end
    # Add initialization of each accumulator
    for a in accus
        push!(functionbody.args,init(a))
    end
    loopbody = Expr(:block,Any[])
    # Add accu step to loop body
    for a in accus
        push!(loopbody.args,accu(a))
    end
    push!(functionbody.args,Expr(:for,:(i=1:n),loopbody))
    # Add returned values
    ret=Expr(:tuple)
    for f in features
        push!(ret.args,final(f))
    end
    push!(functionbody.args,Expr(:return,ret))
    return(Expr(:function,:(getFeatures(e::Extreme,area,lons,lats,prearrays)),functionbody))
end

function getPreArrays(maxlen,ztype,features...)
  arar=Array(Vector{ztype},0)
  for f in features
    if f==Quantile
      l=1
      while l<maxlen
        l=l*10
        push!(arar,Array(ztype,l))
      end
    end
  end
  return(arar)
end


end #Module
