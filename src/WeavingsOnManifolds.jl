module WeavingsOnManifolds
import HomotopyOpt.HomotopyContinuation: Expression, Variable, @var, evaluate, differentiate, System, real_solutions, solve
import HomotopyOpt: ConstraintVariety, findminima
import LinearAlgebra: pinv, norm, cross, nullspace
import Implicit3DPlotting.GLMakiePlottingLibrary: Figure, Axis3, hidespines!, hidedecorations!, Point3f0, scatter!, linesegments!, RGBA
import Implicit3DPlotting: plot_implicit_surface!

export test, computeOptimalWeaving

#TODO angle constraints, repulsive potential between neighboring points (mathematically through paper), part of a minimal surface?, hyperbolic geometry?, PBC?, 6-fold diamond patch, or saddle shape 

struct WeavingOnManifold
    manifold
    offsetList
    offset
    constraints
    coordinateVariables
    bars
    cables
    planes
    variablePositions
    samples
    outsideindices

    function WeavingOnManifold(offsetList::Vector{Bool}, bars::Vector, cables::Vector, planes::Vector; manifold="sphere", implicitManifold=x->x[1]^2+x[2]^2+x[3]^2-1, offset=0.1, torusknot=(true, 3), angleConstraints=true, samples=[])
        (offset<=0.5 && offset>=0.05) || @warn "We recommend an offset between 5% and 50%."
        manifold!="free" || typeof(implicitManifold([1.,2,3]))==Float64

        @var x[1:3, 1:length(offsetList)]

        xs = zeros(Expression,(3,length(offsetList)))
        manifoldEquations, barEquations, angleEquations, planeEquations, xvarz, positions, outsideindices, cableJoints = Vector{Expression}([]), Vector{Expression}([]), Vector{Expression}([]), Vector{Expression}([]), Vector{Variable}([]), [], [], [], []
        
        
        for cable in cables
            if !any(t->(cable[1] in t)||(cable[2] in t), cableJoints)
                push!(cableJoints, [cable[1],cable[2]])
            else
                append!(cableJoints[findfirst(t->(cable[1] in t)||(cable[2] in t), cableJoints)], [cable[1],cable[2]])
            end
        end
        cableJoints = [collect(Set(joints)) for joints in cableJoints]

        for i in findall(t->t,offsetList), j in findall(t->t,offsetList)
            if j<=i
                continue
            end
            push!(outsideindices, (i,j))
        end

        if isequal(lowercase(manifold), "sphere")
            #xs[:,1] = [offsetList[1] ? 1+offset : 1-offset, 0, 0]
            #xs[:,2] = [0, x[2,2], x[3,2]]

            xs[:,1] = isempty(samples) ? (offsetList[1] ? (1+offset) * [ 1,0,0] : (1-offset) * [ 1,0,0]) : samples[1]
            xs[:,2:length(offsetList)] .= x[:,2:length(offsetList)]

            for bar in bars
                xs[:,bar[2]] = offsetList[bar[2]] ? xs[:,bar[1]]*(1+offset)/(1-offset) : xs[:,bar[1]]*(1-offset)/(1+offset)
            end

            for i in 1:size(xs)[2], j in 1:size(xs)[1]
                try
                    push!(xvarz, Variable(xs[j,i]))
                    push!(positions, (j,i))
                catch
                    continue
                end
            end

            samples = isempty(samples) ? append!(samples, [offsetList[1] ? (1+offset) * [ 1,0,0] : (1-offset) * [ 1,0,0]]) : samples
            manifoldEquations = [offsetList[bar[1]] ? sum( xs[:,bar[1]].^2 ) - (1+offset)^2 : sum( xs[:,bar[1]].^2 ) - (1-offset)^2 for bar in bars]
            planeEquations = vcat([[cross(xs[:,plane[2]]-xs[:,plane[1]], xs[:,plane[3]]-xs[:,plane[2]])'*(xs[:,plane[el>length(plane) ? 1 : el]]-xs[:,plane[el-1]]) for el in 4:length(plane)+1] for plane in planes]...)
        elseif isequal(lowercase(manifold), "torus")
            xs[:,1] = [1, 0, offsetList[1] ? 0.5+offset : 0.5+offset]
            #xs[:,2] = [torusknot[1] ? cos(2*pi/torusknot[2])*x[2,2] : 0, torusknot[1] ? sin(2*pi/torusknot[2])*x[2,2] : x[2,2], x[3,2]]
            xs[:,2:length(offsetList)] .= x[:,2:length(offsetList)]
            
            if !torusknot[1]
                append!(barEquations, [sum((xs[:,bar[1]]-xs[:,bar[2]]).^2)-(2*offset)^2 for bar in bars])
            else
                append!(barEquations, [sum((xs[:,bar[1]]-xs[:,bar[2]]).^2)-(1+2*offset*(offsetList[bar[1]]+offsetList[bar[2]]-1))^2 for bar in bars])
            end

            for i in 1:size(xs)[2], j in 1:size(xs)[1]
                try
                    push!(xvarz, Variable(xs[j,i]))
                    push!(positions, (j,i))
                catch
                    continue
                end
            end

            samples = isempty(samples) ? append!(samples, [[1, 0, offsetList[1] ? 0.5+offset : 0.5+offset], [0,0.5,0]]) : samples
            manifoldEquations = [offsetList[i] ? ((1-(0.5+offset)^2)+sum(xs[:,i].^2))^2-4*sum(xs[1:2,i].^2) : ((1-(0.5-offset)^2)+sum(xs[:,i].^2))^2-4*sum(xs[1:2,i].^2) for i in 1:length(offsetList)]
        elseif isequal(lowercase(manifold), "free")            
            xs[:,1] = sols[1]
            xs[:,2:length(offsetList)] .= x[:,2:length(offsetList)]

            append!(barEquations, [sum((xs[:,bar[1]]-xs[:,bar[2]]).^2)-(2*offset)^2 for bar in bars])
            samples = isempty(samples) ? append!(samples, [sols[1],sols[2]]) : samples
            append!(xvarz, x[2:3,2]); append!(positions, [(2,2),(3,2)])
            foreach(i->append!(xvarz, x[:,i]), 3:length(offsetList)); foreach(i->append!(positions, [(1,i),(2,i),(3,i)]), 3:length(offsetList))
            manifoldEquations = [implicitManifold(xs[:,i])+offset*2*(offsetList[i]-0.5) for i in 2:length(offsetList)]
        else
            throw(error("Only weavings on the sphere are implemented right now!"))
        end
        
        totalContactCombinatorics = [[map(t->t[1]==bar[1] ? t : (t[2],t[1]), cables[findall(cable->bar[1] in cable, cables)]), bar, map(t->t[1]==bar[1] ? t : (t[2],t[1]), cables[findall(cable->bar[2] in cable, cables)])] for bar in bars]
        angleEquations = angleConstraints ? vcat([[((xs[:,contact[1][1][1]]-xs[:,contact[1][1][2]])'*(xs[:,contact[2][1]]-xs[:,contact[2][2]]))^2*((xs[:,contact[1][2][2]]-xs[:,contact[1][2][1]])'*(xs[:,contact[1][2][2]]-xs[:,contact[1][2][1]]))-((xs[:,contact[1][2][2]]-xs[:,contact[1][2][1]])'*(xs[:,contact[2][1]]-xs[:,contact[2][2]]))^2*((xs[:,contact[1][1][1]]-xs[:,contact[1][1][2]])'*(xs[:,contact[1][1][1]]-xs[:,contact[1][1][2]])), ((xs[:,contact[3][1][1]]-xs[:,contact[3][1][2]])'*(xs[:,contact[2][1]]-xs[:,contact[2][2]]))^2*((xs[:,contact[3][2][2]]-xs[:,contact[3][2][1]])'*(xs[:,contact[3][2][2]]-xs[:,contact[3][2][1]]))-((xs[:,contact[3][2][2]]-xs[:,contact[3][2][1]])'*(xs[:,contact[2][1]]-xs[:,contact[2][2]]))^2*((xs[:,contact[3][1][1]]-xs[:,contact[3][1][2]])'*(xs[:,contact[3][1][1]]-xs[:,contact[3][1][2]]))] for contact in totalContactCombinatorics]...) : angleEquations

        energyFunction = sum([sum((xs[:,cable[1]] - xs[:,cable[2]]).^2) for cable in cables])
        new(lowercase(manifold), offsetList, offset, Vector{Expression}(vcat(manifoldEquations, barEquations, planeEquations, angleEquations)), xvarz, bars, cables, planes, positions, samples, outsideindices)
    end
end

function toMatrix(configuration, Weave::WeavingOnManifold)
    p0 = zeros(Number,3,Int(length(Weave.offsetList)))
    global count = 1
    p0[:,1] = Weave.samples[1]

    for pos in Weave.variablePositions
        p0[pos[1],pos[2]] = configuration[count]
        global count = count+1
    end

    if Weave.manifold=="sphere"
        for bar in Weave.bars
            p0[:,bar[2]] = Weave.offsetList[bar[2]] ? p0[:,bar[1]]*(1+Weave.offset)/(1-Weave.offset) : p0[:,bar[1]]*(1-Weave.offset)/(1+Weave.offset)
        end
    end
    return(p0)
end

function toArray(p, Weave::WeavingOnManifold)
    configuration = Vector{Float64}([])
    positions = Weave.variablePositions
    foreach(pos->append!(configuration, p[pos[1], pos[2]]), positions)
    return configuration
end

function energyFunction(configuration, Weave::WeavingOnManifold)
    p = toMatrix(configuration, Weave)
    Q = 2*sum([ sum( (p[:,cable[1]] - p[:,cable[2]]).^2 ) for cable in Weave.cables])
    Q += -sum([ sum( (p[:,t[1]] - p[:,t[2]]).^2 ) for t in Weave.outsideindices])
    return Q
end

function plotWeaving(configuration::Vector{Float64}, Weave::WeavingOnManifold)
    p0 = toMatrix(configuration, Weave)
    fig = Figure(size = (1000,1000))
    ax = Axis3(fig[1,1], aspect=(1.,1,1))
   # hidespines!(ax); hidedecorations!(ax);
    implicit_manifold = x->(Weave.manifold == "torus") ? (0.75 + (x[1]^2+x[2]^2+x[3]^2))^2 - 4*(x[1]^2+x[2]^2) : x[1]^2+x[2]^2+x[3]^2-1
    plot_implicit_surface!(ax, implicit_manifold; wireframe=false, transparency=true, color=RGBA(0.75,0.75,0.75,0.6), samples = (75,75,75))
    foreach(bar->linesegments!(ax, [Point3f0(p0[:,bar[1]]), Point3f0(p0[:,bar[2]])]; color=:blue, linewidth=8), Weave.bars)
    foreach(cable->linesegments!(ax, [Point3f0(p0[:,cable[1]]), Point3f0(p0[:,cable[2]])]; color=:red, linewidth=8), Weave.cables)

    scatter!(ax, [Point3f0(p0[:,i]) for i in 1:size(p0)[2]]; color=:black, markersize=35)
    display(fig)
end

function test_sphere_a()
    Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true, true,false,true,false, false,true,false,true], [(1,5),(2,11),(3,7),(4,9),(6,10),(8,12)], [(1,2),(2,3),(3,4),(4,1), (5,6),(6,7),(7,8),(8,5), (9,10),(10,11),(11,12),(12,9)], [(1,2,3,4), (5,6,7,8), (9,10,11,12)])
    p0 = [1 0 0; 0 -1 0; -1 0 0; 0 1 0; 1 0 0; 0 0 -1; -1 0 0; 0 0 1; 0 1 0; 0 0 -1; 0 -1 0; 0 0 1]'
    initialConfiguration = toArray(p0, Weave) #+ (randn(Float64, length(toArray(p0, Weave))) .- 0.5)*0.05
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave)
    q = computeOptimalWeaving(q, Weave)
    normalspace=evaluate.(differentiate(Weave.constraints, Weave.coordinateVariables), Weave.coordinateVariables=>q)
    display("Tangent direction")
    display(toMatrix(collect(nullspace(normalspace)[1:end,1]), Weave))
    plotWeaving(q, Weave)
end

function test_sphere_a2()
    Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true,false,true,false,true, true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true, true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true], 
            [(4, 43), (9, 28), (15, 47), (29, 38), (8, 40), (2, 19), (13, 23), (12, 44), (21, 37), (26, 33), (30, 46), (22, 45), (11, 36), (7, 48), (3, 35), (14, 31), (17, 42), (6, 32), (16, 39), (5, 24), (18, 34), (10, 20), (25, 41), (1, 27)], 
            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,1), (9,10),(10,11),(11,12),(12,13),(13,14),(14,15),(15,16),(16,9), (17,18),(18,19),(19,20),(20,21),(21,22),(22,23),(23,24),(24,17), (25,26),(26,27),(27,28),(28,29),(29,30),(30,31),(31,32),(32,25), (33,34),(34,35),(35,36),(36,37),(37,38),(38,39),(39,40),(40,33), (41,42),(42,43),(43,44),(44,45),(45,46),(46,47),(47,48),(48,41)], []; offset = 0.1, samples = [0.9*[-1/3, sqrt(7/9), 1/3], 1.1*[1/3, sqrt(7/9), 1/3]])
    p0 = [-1/3 sqrt(7/9) 1/3; 1/3 sqrt(7/9) 1/3; sqrt(7/9) 1/3 1/3; sqrt(7/9) -1/3 1/3; 1/3 -sqrt(7/9) 1/3; -1/3 -sqrt(7/9) 1/3; -sqrt(7/9) -1/3 1/3; -sqrt(7/9) 1/3 1/3;
    -1/3 sqrt(7/9) -1/3; 1/3 sqrt(7/9) -1/3; sqrt(7/9) 1/3 -1/3; sqrt(7/9) -1/3 -1/3; 1/3 -sqrt(7/9) -1/3; -1/3 -sqrt(7/9) -1/3; -sqrt(7/9) -1/3 -1/3; -sqrt(7/9) 1/3 -1/3;
    1/3 -1/3 sqrt(7/9); 1/3 1/3 sqrt(7/9); 1/3 sqrt(7/9) 1/3; 1/3 sqrt(7/9) -1/3; 1/3 1/3 -sqrt(7/9); 1/3 -1/3 -sqrt(7/9); 1/3 -sqrt(7/9) -1/3; 1/3 -sqrt(7/9) 1/3;
    -1/3 -1/3 sqrt(7/9); -1/3 1/3 sqrt(7/9); -1/3 sqrt(7/9) 1/3; -1/3 sqrt(7/9) -1/3; -1/3 1/3 -sqrt(7/9); -1/3 -1/3 -sqrt(7/9); -1/3 -sqrt(7/9) -1/3; -1/3 -sqrt(7/9) 1/3;
    -1/3 1/3 sqrt(7/9); 1/3 1/3 sqrt(7/9); sqrt(7/9) 1/3 1/3; sqrt(7/9) 1/3 -1/3; 1/3 1/3 -sqrt(7/9); -1/3 1/3 -sqrt(7/9); -sqrt(7/9) 1/3 -1/3; -sqrt(7/9) 1/3 1/3;
    -1/3 -1/3 sqrt(7/9); 1/3 -1/3 sqrt(7/9); sqrt(7/9) -1/3 1/3; sqrt(7/9) -1/3 -1/3; 1/3 -1/3 -sqrt(7/9); -1/3 -1/3 -sqrt(7/9); -sqrt(7/9) -1/3 -1/3; -sqrt(7/9) -1/3 1/3;
    ]'

    #=
    bars = []  
    for i in 1:size(p0)[2], j in i+1:size(p0)[2]
        if norm(p0[:,i]-p0[:,j])<=10^(-10)
            push!(bars,(i,j))
        end
    end 
    println(collect(Set(bars)))=#

    initialConfiguration = toArray(p0, Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving( q, Weave )
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave)
end

function test_sphere_b()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false, false,true,false,true,false,true, false,true,false,true,false,true, true,false,true,false,true,false], [(1,7), (4,10), (14,20), (17,23), (3,15), (6,18), (9,16), (12,13), (11,19), (5,24), (2,21), (8,22)], [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1), (7,8),(8,9),(9,10),(10,11),(11,12),(12,7), (13,14),(14,15),(15,16),(16,17),(17,18),(18,13), (19,20),(20,21),(21,22),(22,23),(23,24),(24,19)], [(1,2,3,4,5,6), (7,8,9,10,11,12), (13,14,15,16,17,18), (19,20,21,22,23,24)])
    p0 = [1 0 0; cos(2*pi/6) sin(2*pi/6) 0; cos(4*pi/6) sin(4*pi/6) 0; cos(6*pi/6) sin(6*pi/6) 0; cos(8*pi/6) sin(8*pi/6) 0; cos(10*pi/6) sin(10*pi/6) 0;
        1 0 0; cos(2*pi/6) 0 sin(2*pi/6); cos(4*pi/6) 0 sin(4*pi/6); cos(6*pi/6) 0 sin(6*pi/6); cos(8*pi/6) 0 sin(8*pi/6); cos(10*pi/6) 0 sin(10*pi/6);
        cos(10*pi/6) 0 sin(10*pi/6); 0 sqrt(1/2) -sqrt(1/2); cos(4*pi/6) sin(4*pi/6) 0; cos(4*pi/6) 0 sin(4*pi/6); 0 -sqrt(1/2) sqrt(1/2); cos(10*pi/6) sin(10*pi/6) 0;
        cos(8*pi/6) 0 sin(8*pi/6); 0 sqrt(1/2) -sqrt(1/2); cos(2*pi/6) sin(2*pi/6) 0; cos(2*pi/6) 0 sin(2*pi/6); 0 -sqrt(1/2) sqrt(1/2); cos(8*pi/6) sin(8*pi/6) 0]'
    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave)
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave)
end

function test_torus_trefoil()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true, true, true, true, true, true], [(1,4), (2,5), (3,6)], [(1,2), (2,3), (3,4), (4,5), (5,6), (6,1)], []; manifold="torus", offset=0.)
    p0 = [1 0 0.5; cos(2*pi/3) sin(2*pi/3) -0.5; cos(4*pi/3) sin(4*pi/3) 0.5; 1 0 -0.5; cos(2*pi/3) sin(2*pi/3) 0.5; cos(4*pi/3) sin(4*pi/3) -0.5]'
    
    initialConfiguration = toArray(p0, Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave)
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave)
end


function newtonCorrect(q::Vector{Float64}, variables, equations; tol = 1e-8)
    global qnew, damping = q, 0.25
	jac = hcat([differentiate(eq, variables) for eq in equations]...)
    while( norm(evaluate.(equations, variables=>q)) > tol )
		J = Matrix{Float64}(evaluate.(jac, variables=>q))
		global qnew = q .- damping*pinv(J)'*evaluate.(equations, variables=>q)
        display(norm(evaluate.(equations, variables=>q)))
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
    resultmin = findminima(q, 1e-3, G, Q; whichstep="gaussnewtonstep", maxseconds=100, stepdirection = "gradientdescent")
    return(resultmin.computedpoints[end])
end

test_torus_trefoil()
#test_torus_trefoil()
end
