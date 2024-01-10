module WeavingsOnManifolds
import HomotopyContinuation: Expression, Variable, @var
import Implicit3DPlotting.GLMakiePlottingLibrary: Figure, Axis3, hidespines!, hidedecorations!, Point3f0, scatter!, linesegments!, RGBA
import Implicit3DPlotting: plot_implicit_surface!

mutable struct WeavingOnManifold
    manifold
    constraints
    energyFunction
    lagrangeSystem
    coordinateVariables
    lagrangeMultipliers
    bars
    cables

    function WeavingOnManifold(offsetList::Vector{Bool}, bars::Vector, cables::Vector, manifold::String, offset::Float64)
        (offset<=0.5 && offset>=0.05) || @warn "We recommend an offset between 5% and 50%."
        @var xs[1:3, 1:length(offsetList)]
        manifoldEquations, barEquations = [], []
        if isequal(lowercase(manifold), "sphere")
            manifoldEquations = [offsetList[i] ? sum( xs[:,i].^2 ) - (1+offset)^2 : sum( xs[:,i].^2 ) - (1-offset)^2 for i in 1:length(offsetList)]
        else
            throw(error("Only weavings on the sphere are implemented right now!"))
        end
        barEquations = [sum( (xs[:,bar[1]] - xs[:,bar[2]]).^2 ) - (2*offset)^2 for bar in bars]
        energyFunction = sum([sum((xs[:,cable[1]] - xs[:,cable[2]]).^2) for cable in cables])
        @var λ[1:(length(manifoldEquations)+length(barEquations))]
        new(lowercase(manifold), vcat(manifoldEquations, barEquations), energyFunction, energyFunction + λ'*vcat(manifoldEquations, barEquations), reduce(vcat, xs), λ, bars, cables)
    end
end

function toMatrix(configuration::Vector{Float64})
    p0 = zeros(Float64,3,Int(length(configuration)/3))
    global count = 1
    for i in 1:Int(length(configuration)/3)
        p0[1:3,i] = configuration[count:count+2]
        global count = count+3
    end
    return(p0)
end

function plotWeaving(configuration::Vector{Float64}, Weave::WeavingOnManifold)
    p0 = toMatrix(configuration)
    fig = Figure(size = (1000,1000))
    ax = Axis3(fig[1,1], aspect=(1.,1,1))
    hidespines!(ax)
    hidedecorations!(ax)
    plot_implicit_surface!(ax, x->x[1]^2+x[2]^2+x[3]^2-1; wireframe=false, transparency=true, color=RGBA(0.5,0.5,0.5,0.6))
    foreach(bar->linesegments!(ax, [Point3f0(p0[:,bar[1]]), Point3f0(p0[:,bar[2]])]; color=:blue, linewidth=8), Weave.bars)
    foreach(cable->linesegments!(ax, [Point3f0(p0[:,cable[1]]), Point3f0(p0[:,cable[2]])]; color=:red, linewidth=8), Weave.cables)

    scatter!(ax, [Point3f0(p0[:,i]) for i in 1:size(p0)[2]]; color=:black, markersize=35)
    display(fig)
end

function test()
    Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true],[],[(1,2),(2,3),(3,4),(4,1)],"sphere",0.1)
    plotWeaving([1+0.1,0,0,0,-1+0.1,0,-1-0.1,0,0,0,1-0.1,0], Weave::WeavingOnManifold)
end

function computeOptimalWeaving( Weaving::WeavingOnManifold)

end

end
