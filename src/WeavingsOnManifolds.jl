module WeavingsOnManifolds
import HomotopyContinuation: Expression, Variable, @var, evaluate, differentiate
import LinearAlgebra: pinv, norm, cross
import Implicit3DPlotting.GLMakiePlottingLibrary: Figure, Axis3, hidespines!, hidedecorations!, Point3f0, scatter!, linesegments!, RGBA
import Implicit3DPlotting: plot_implicit_surface!
import HomotopyOpt: ConstraintVariety, findminima

export test, computeOptimalWeaving

struct WeavingOnManifold
    manifold
    offsetList
    offset
    constraints
    coordinateVariables
    bars
    cables
    planes

    function WeavingOnManifold(offsetList::Vector{Bool}, bars::Vector, cables::Vector, planes::Vector, manifold::String, offset::Float64)
        (offset<=0.5 && offset>=0.05) || @warn "We recommend an offset between 5% and 50%."
        @var x[1:3, 1:length(offsetList)]

        xs = zeros(Expression,(3,length(offsetList)))
        manifoldEquations, barEquations = [], []
        if isequal(lowercase(manifold), "sphere")
            xs[:,1] = [offsetList[1] ? 1+offset : 1-offset, 0, 0]
            xs[:,2] = [0, x[2,2], x[3,2]]
            for i in 3:length(offsetList)
                xs[:,i] = x[:,i]
            end

            for bar in bars
                xs[:,bar[2]] = offsetList[bar[2]] ? xs[:,bar[1]]*(1+offset)/(1-offset) : xs[:,bar[1]]*(1-offset)/(1+offset)
            end
            for i in 1:size(xs)[2]
                for j in 1:size(xs)[1]
            manifoldEquations = [offsetList[i] ? sum( xs[:,i].^2 ) - (1+offset)^2 : sum( xs[:,i].^2 ) - (1-offset)^2 for i in 1:length(offsetList)]
        else
            throw(error("Only weavings on the sphere are implemented right now!"))
        end
        #barEquations = [sum( (xs[:,bar[1]] - xs[:,bar[2]]).^2 ) - (2*offset)^2 for bar in bars]
        planeEquations = vcat([[cross(xs[:,plane[2]]-xs[:,plane[1]], xs[:,plane[3]]-xs[:,plane[2]])'*(xs[:,plane[4]]-xs[:,plane[3]]), cross(xs[:,plane[2]]-xs[:,plane[1]], xs[:,plane[3]]-xs[:,plane[2]])'*(xs[:,plane[1]]-xs[:,plane[4]])] for plane in planes]...)
        energyFunction = sum([sum((xs[:,cable[1]] - xs[:,cable[2]]).^2) for cable in cables])
        new(lowercase(manifold), offsetList, offset, Vector{Expression}(vcat(manifoldEquations, barEquations, planeEquations)), #=insert x coordinates here=#, bars, cables, planes)
    end
end

function toMatrix(configuration, Weave::WeavingOnManifold)
    p0 = zeros(Number,3,Int((length(configuration)+4)/3))
    if isequal(Weave.manifold, "sphere")
        global count = 3
        p0[:,1] = [Weave.offsetList[1] ? 1+Weave.offset : 1-Weave.offset, 0, 0]
        p0[:,2] = [0, configuration[1], configuration[2]]
    else
        throw(error("Only weavings on the sphere are implemented right now!"))
    end

    for i in 3:size(p0)[2]
        p0[1:3,i] = configuration[count:count+2]
        global count = count+3
    end
    return(p0)
end

function toArray(p)
    configuration = Vector{Float64}([])
    append!(configuration, p[2:3,2])
    foreach(i->append!(configuration, p[:,i]), 3:size(p)[2])
    return configuration
end

function energyFunction(configuration, Weave::WeavingOnManifold)
    Q = 0
    p = toMatrix(configuration, Weave)
    for cable in Weave.cables
        Q += sum( (p[:,cable[1]] - p[:,cable[2]]).^2 )
    end
    return Q
end

function plotWeaving(configuration::Vector{Float64}, Weave::WeavingOnManifold)
    p0 = toMatrix(configuration, Weave)
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
    Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true, true,false,true,false, false,true,false,true], [(1,5),(2,11),(3,7),(4,9),(6,10),(8,12)], [(1,2),(2,3),(3,4),(4,1), (5,6),(6,7),(7,8),(8,5), (9,10),(10,11),(11,12),(12,9)], [(1,2,3,4), (5,6,7,8), (9,10,11,12)], "sphere",0.1)
    p0 = [1 0 0; 0 -1 0; -1 0 0; 0 1 0; 1 0 0; 0 0 -1; -1 0 0; 0 0 1; 0 1 0; 0 0 -1; 0 -1 0; 0 0 1]'
    display(Weave)
    initialConfiguration = toArray(p0)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave)
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave)
end

function newtonCorrect(q::Vector{Float64}, variables, equations; tol = 1e-8)
    global damping = 0.5
	global qnew = q
	jac = Base.hcat([differentiate(eq, variables) for eq in equations]...)
    while(norm(evaluate.(equations, variables=>q)) > tol)
		J = Matrix{Float64}(evaluate.(jac, variables=>q))
		global qnew = q .- damping*pinv(J)'*evaluate.(equations, variables=>q)
		if norm(evaluate.(equations, variables=>qnew)) <= norm(evaluate.(equations, variables=>q))
			global damping = damping*1.2
		else
			global damping = damping/2
		end
		q = qnew
	end
	return q
end

function computeOptimalWeaving(initialConfiguration::Vector{Float64}, Weave::WeavingOnManifold; tol = 1e-4)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints)
    Q = x->energyFunction(x, Weave)
    G = ConstraintVariety(Weave.coordinateVariables, Weave.constraints, length(Weave.coordinateVariables), length(Weave.coordinateVariables)-length(Weave.constraints))
    resultmin = findminima(q, 1e-3, G, Q; whichstep="gaussnewtonstep", maxseconds=25)
    return(resultmin.computedpoints[end])
end

test()

end
