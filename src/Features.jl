module Features

#Define Accumulator behaviors
abstract Accumulator
init(t::Type{Accumulator}) = error("Init behavior for Accumulator $t not defined.")
accu(t::Type{Accumulator}) = error("Accumulation behavior for Accumulator $t not defined.")
type Max_z_acc <: Accumulator end
init(::Type{Max_z_acc})= :(max_z_acc=typemin(eltype(z)))
accu(::Type{Max_z_acc})= :(if z[i]>max_z_acc max_z_acc=z[i] end)
type Min_z_acc <: Accumulator end
init(::Type{Max_z_acc})= :(min_z_acc=typemax(eltype(z)))
accu(::Type{Max_z_acc})= :(if z[i]<min_z_acc min_z_acc=z[i] end)
type Sum_z_acc <: Accumulator end
init(::Type{Sum_z_acc})= :(sum_z_acc=zero(eltype(z)))
accu(::Type{Sum_z_acc})= :(sum_z_acc=sum_z_acc+z[i])

type Mean end
accumulators(::Type{Mean})=(Sum_z_acc,)
final(::Type{Mean})= :(sum_z_acc/length(z))

function calcFeatureFunction(features::Type...)
    # Get required accumulators
    accus = Array(Any,0)
    for feat in features
        for a in accumulators(feat)
            in(a,accus) || push!(accus,a)
        end
    end
    functionbody = quote
        x=e.locs
        z=e.zvalues
        n=length(z)
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
    return(Expr(:function,:(getFeatures(e::Extreme,)),functionbody))
end



end #Module
