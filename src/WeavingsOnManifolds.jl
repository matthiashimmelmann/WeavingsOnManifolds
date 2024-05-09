module WeavingsOnManifolds
import HomotopyOpt.HomotopyContinuation: Expression, Variable, @var, evaluate, differentiate, System, real_solutions, solve
import HomotopyOpt: ConstraintVariety, findminima
import LinearAlgebra: pinv, norm, cross, nullspace
import Implicit3DPlotting.GLMakiePlottingLibrary: hidedecorations!, hidespines!, Figure, Axis3, hidespines!, hidedecorations!, Point3f0, scatter!, linesegments!, RGBA
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
    periodicBoundary

    function WeavingOnManifold(offsetList::Vector{Bool}, bars::Vector, cables::Vector, planes::Vector; manifold="sphere", implicitManifold=x->x[1]^2+x[2]^2+x[3]^2-1, offset=0.1, torusknot=(true, 3), angleConstraints=false, samples=[], PBC=[])
        (offset<=0.5 && offset>=0.05) || @warn "We recommend an offset between 5% and 50%."
        manifold!="free" || typeof(implicitManifold([1.,2,3]))==Float64
        ([i for i in 1:length(offsetList)] == sort(vcat([bar[1] for bar in bars], [bar[2] for bar in bars])) && (Set(1:length(offsetList)) == Set(vcat([cable[1] for cable in cables], [cable[2] for cable in cables])))) || throw(error("Cables and Bars don't cover every vertex"))

        @var x[1:3, 1:length(offsetList)]

        xs = zeros(Expression,(3,length(offsetList)))
        manifoldEquations, barEquations, angleEquations, planeEquations, periodicBoundaryConditions, xvarz, positions, outsideindices, cableJoints = Vector{Expression}([]), Vector{Expression}([]), Vector{Expression}([]), Vector{Expression}([]), Vector{Expression}([]), Vector{Variable}([]), [], [], []
        
        
        
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
            planeEquations = vcat([[cross(xs[:,plane[2]]-xs[:,plane[1]], xs[:,plane[4]]-xs[:,plane[3]])'*(xs[:,plane[el>length(plane) ? 1 : el]]-xs[:,plane[el-1]]) for el in vcat(3,5:length(plane)+1)] for plane in planes]...)
        elseif isequal(lowercase(manifold), "flattorus")
            xs[:,1] = samples[1]+[0,0,offsetList[1] ? offset : -offset]
            xs[:,2:length(offsetList)] .= x[:,2:length(offsetList)]
            bars = map(bar->offsetList[bar[1]] ? bar : (bar[2],bar[1]), bars)

            for i in 1:size(xs)[2], j in 1:size(xs)[1]
                push!(xvarz, x[j,i])
                push!(positions, (j,i))
            end

            outsideindices = [(1,5), (5,9), (5,7), (3,5), (13,11),(13,17),(11,15),(15,17)]

            for cond in PBC
                if offsetList[cond[1]]    
                    xs[Float64(cond[3][1])==1. ? 1 : 2, cond[1]] = 1.
                    index1 = findfirst(pos->pos[2]==cond[1]&&((Float64(cond[3][1])==1. && pos[1]==1)||(Float64(cond[3][2])==1. && pos[1]==2)), positions)
                    deleteat!(positions,index1); deleteat!(xvarz,index1)

                    xs[:, cond[2]] = xs[:, cond[1]]-cond[3]
                    index2 = findfirst(pos->pos[2]==cond[2]&&pos[1]==1, positions)
                    if index2!=nothing
                        deleteat!(positions,index2); deleteat!(xvarz,index2)
                    end
                    
                    index3 = findfirst(pos->pos[2]==cond[2]&&pos[1]==2, positions)
                    if index3!=nothing
                        deleteat!(positions,index3); deleteat!(xvarz,index3)
                    end
                    index4 = findfirst(pos->pos[2]==cond[2]&&pos[1]==3, positions)
                    if index4!=nothing
                        deleteat!(positions,index4); deleteat!(xvarz,index4)
                    end 
                else
                    top_bar1 = bars[findfirst(t->cond[1] in t, bars)][1]
                    println(cond[1],", ",top_bar1)
                    xs[Float64(cond[3][1])==1. ? 1 : 2, top_bar1] = 1.
                    index1 = findfirst(pos->pos[2]==top_bar1&&((Float64(cond[3][1])==1. && pos[1]==1)||(Float64(cond[3][2])==1. && pos[1]==2)), positions)
                    if index1!=nothing
                        deleteat!(positions,index1); deleteat!(xvarz,index1)
                    end

                    top_bar2 = bars[findfirst(t->cond[2] in t, bars)][1]
                    xs[:, top_bar2] = xs[:, top_bar1]-cond[3]
                    index2 = findfirst(pos->pos[2]==top_bar2&&pos[1]==1, positions)
                    if index2!=nothing
                        deleteat!(positions,index2); deleteat!(xvarz,index2)
                    end
                    
                    index3 = findfirst(pos->pos[2]==top_bar2&&pos[1]==2, positions)
                    if index3!=nothing
                        deleteat!(positions,index3); deleteat!(xvarz,index3)
                    end
                    index4 = findfirst(pos->pos[2]==top_bar2&&pos[1]==3, positions)
                    if index4!=nothing
                        deleteat!(positions,index4); deleteat!(xvarz,index4)
                    end 
                end
            end

            for bar in bars
                index1 = findfirst(pos->pos[2]==bar[2]&&pos[1]==1, positions)

                if index1!=nothing
                    deleteat!(positions,index1); deleteat!(xvarz,index1)
                    xs[1,bar[2]] = xs[1,bar[1]]
                end

                index2 = findfirst(pos->pos[2]==bar[2]&&pos[1]==2, positions)
                if index2!=nothing
                    deleteat!(positions,index2); deleteat!(xvarz,index2)
                    xs[2,bar[2]] = xs[2,bar[1]]
                end

                index3 = findfirst(pos->pos[2]==bar[2]&&pos[1]==3, positions)
                if index3!=nothing
                    deleteat!(positions,index3); deleteat!(xvarz,index3)
                    xs[3,bar[2]] = xs[3,bar[1]]-2*offset
                end
            end

            #=for cond in PBC
                append!(periodicBoundaryConditions, Vector{Expression}(xs[:,cond[1]]-xs[:,cond[2]]-cond[3]))
            end=#

            samples = [samples[1]+[0,0, offsetList[1] ? offset : -offset]]
            manifoldEquations = [offsetList[i] ? xs[3,i]-offset : xs[3,i]+offset for i in 1:length(offsetList)]
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
        #angleEquations = angleConstraints ? vcat([[((xs[:,contact[1][1][1]]-xs[:,contact[1][1][2]])'*(xs[:,contact[2][1]]-xs[:,contact[2][2]]))^2*((xs[:,contact[1][2][2]]-xs[:,contact[1][2][1]])'*(xs[:,contact[1][2][2]]-xs[:,contact[1][2][1]]))-((xs[:,contact[1][2][2]]-xs[:,contact[1][2][1]])'*(xs[:,contact[2][1]]-xs[:,contact[2][2]]))^2*((xs[:,contact[1][1][1]]-xs[:,contact[1][1][2]])'*(xs[:,contact[1][1][1]]-xs[:,contact[1][1][2]])), ((xs[:,contact[3][1][1]]-xs[:,contact[3][1][2]])'*(xs[:,contact[2][1]]-xs[:,contact[2][2]]))^2*((xs[:,contact[3][2][2]]-xs[:,contact[3][2][1]])'*(xs[:,contact[3][2][2]]-xs[:,contact[3][2][1]]))-((xs[:,contact[3][2][2]]-xs[:,contact[3][2][1]])'*(xs[:,contact[2][1]]-xs[:,contact[2][2]]))^2*((xs[:,contact[3][1][1]]-xs[:,contact[3][1][2]])'*(xs[:,contact[3][1][1]]-xs[:,contact[3][1][2]]))] for contact in totalContactCombinatorics]...) : angleEquations

        energyFunction = sum([sum((xs[:,cable[1]] - xs[:,cable[2]]).^2) for cable in cables])
        new(lowercase(manifold), offsetList, offset, filter(t->t!=0,Vector{Expression}(vcat(manifoldEquations, barEquations, planeEquations, angleEquations, periodicBoundaryConditions))), Vector{Variable}(xvarz), bars, cables, planes, positions, samples, outsideindices, PBC)
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
    elseif Weave.manifold=="flattorus"
        for cond in Weave.periodicBoundary
            if Weave.offsetList[cond[1]]
                p0[Float64(cond[3][1])==1. ? 1 : 2, cond[1]] = 1.
                p0[:, cond[2]] = p0[:, cond[1]]-cond[3]
            else
                top_bar1 = Weave.bars[findfirst(t->cond[1] in t, Weave.bars)][1]
                p0[Float64(cond[3][1])==1. ? 1 : 2, top_bar1] = 1.
                top_bar2 = Weave.bars[findfirst(t->cond[2] in t, Weave.bars)][1]
                p0[:, top_bar2] = p0[:, top_bar1]-cond[3]
            end            
        end

        for bar in Weave.bars
            p0[:,bar[2]] = p0[:,bar[1]]-[0,0,2*Weave.offset]
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
    Q = sum([ sum( (p[:,cable[1]] - p[:,cable[2]]).^2 ) for cable in Weave.cables])
    Q += sum([ 1/sum( (p[:,t[1]] - p[:,t[2]]).^2 ) for t in Weave.outsideindices])
    return Q
end

function plotWeaving(configuration::Vector{Float64}, Weave::WeavingOnManifold; colorscheme=[])
    p0 = toMatrix(configuration, Weave)
    fig = Figure(size = (1000,1000))
    ax = Axis3(fig[1,1], aspect=(1.,1,1))
    hidedecorations!(ax)
    hidespines!(ax)
   # hidespines!(ax); hidedecorations!(ax);
    implicit_manifold = x->(Weave.manifold == "torus") ? (0.75 + (x[1]^2+x[2]^2+x[3]^2))^2 - 4*(x[1]^2+x[2]^2) : ((Weave.manifold == "flattorus") ? x[3] : x[1]^2+x[2]^2+x[3]^2-1)
    plot_implicit_surface!(ax, implicit_manifold; wireframe=false, transparency=true, color=RGBA(0.75,0.75,0.75,0.6), samples = (75,75,75), xlims = (Weave.manifold == "flattorus") ? (0,1) : (-1-Weave.offset, 1+Weave.offset), ylims = (Weave.manifold == "flattorus") ? (0,1) : (-1-Weave.offset, 1+Weave.offset), zlims = (Weave.manifold == "flattorus") ? (-2*Weave.offset-0.01,2*Weave.offset+0.01) : (-1-Weave.offset, 1+Weave.offset))
    foreach(bar->linesegments!(ax, [Point3f0(p0[:,bar[1]]), Point3f0(p0[:,bar[2]])]; color=:blue, linewidth=8), Weave.bars)

    if isempty(colorscheme)
        foreach(cable->linesegments!(ax, [Point3f0(p0[:,cable[1]]), Point3f0(p0[:,cable[2]])]; color=:red, linewidth=8), Weave.cables)
    else
        cableJoints = []
        for cable in Weave.cables
            if !any(t->(cable[1] in t)||(cable[2] in t), cableJoints)
                push!(cableJoints, [cable[1],cable[2]])
            else
                append!(cableJoints[findfirst(t->(cable[1] in t)||(cable[2] in t), cableJoints)], [cable[1],cable[2]])
            end
        end
        foreach(cable->linesegments!(ax, [Point3f0(p0[:,cable[1]]), Point3f0(p0[:,cable[2]])]; color=colorscheme[findfirst(t->cable[1] in t, cableJoints)], linewidth=8), Weave.cables)
    end

    scatter!(ax, [Point3f0(p0[:,i]) for i in 1:size(p0)[2]]; color=:black, markersize=35)
    display(fig)
end

function toSphere(configuration, Weave)
    p0 = toMatrix(configuration, Weave)
    newEdgeList = []
    newPointList = []
    display(p0)
    for i in 1:size(p0)[2]
        #p0[:,i] = p0[:,i]./norm(p0[:,i])
    end
    display(p0)
    outer = 1
    epsilon = 0.15
    for edge in Weave.cables
        outer = norm(p0[:,edge[1]])>1 ? 1 : -1
        for t in 0:0.025:0.975
            curP = t*(p0[:,edge[2]]-p0[:,edge[1]])+p0[:,edge[1]]
            #push!(newPointList, (outer==1 ? (1+epsilon-(6*epsilon)*t^2+(4*epsilon)*t^3) : (1+epsilon-(6*epsilon)*(1-t)^2+(4*epsilon)*(1-t)^3))*curP./norm(curP))
            push!(newPointList, curP ./ norm(curP))
            push!(newEdgeList, (length(newPointList), length(newPointList)+1))
        end
        curP = p0[:,edge[2]]
        push!(newPointList, #=(1-outer*epsilon)*=#curP ./ norm(curP))
    end
    for i in size(p0)[2]
        push!(newPointList, p0[:,i] ./ norm(p0[:,i]))
    end
    open("weavingframeworkmodelonsphereb2.poly", "w") do f
        write(f, "POINTS\n")
        foreach(i->write(f, string("$(i): ", newPointList[i][1], " ", newPointList[i][2], " ", newPointList[i][3],"\n")), 1:length(newPointList))
        write(f,"POLYS\n")
        foreach(i->write(f, string("$(i): ", newEdgeList[i][1], " ", newEdgeList[i][2],"\n")), 1:length(newEdgeList))
        write(f,"END")
    end
end

function test_sphere_a()
    Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true, true,false,true,false, false,true,false,true], [(1,5),(2,11),(3,7),(4,9),(6,10),(8,12)], [(1,2),(2,3),(3,4),(4,1), (5,6),(6,7),(7,8),(8,5), (9,10),(10,11),(11,12),(12,9)], []#=[(1,2,3,4), (5,6,7,8), (9,10,11,12)]=#)
    p0 = [1 0 0; 0 -1 0; -1 0 0; 0 1 0; 1 0 0; 0 0 -1; -1 0 0; 0 0 1; 0 1 0; 0 0 -1; 0 -1 0; 0 0 1]'
    initialConfiguration = toArray(p0, Weave) #+ (randn(Float64, length(toArray(p0, Weave))) .- 0.5)*0.05
    display(Weave.coordinateVariables)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-12)
    q = computeOptimalWeaving(q, Weave)
    normalspace=evaluate.(differentiate(Weave.constraints, Weave.coordinateVariables), Weave.coordinateVariables=>q)
    display("Tangent direction")
    display(toMatrix(collect(nullspace(normalspace)[1:end,1]), Weave))
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:green,:blue,:red])
    toSphere(q, Weave)
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
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave; colorscheme=[:yellow,:yellow,:cyan,:cyan,:purple,:purple])
    toSphere(q, Weave)
end

function test_sphere_b()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false, false,true,false,true,false,true, false,true,false,true,false,true, true,false,true,false,true,false], [(1,7), (4,10), (14,20), (17,23), (3,15), (6,18), (9,16), (12,13), (11,19), (5,24), (2,21), (8,22)], [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1), (7,8),(8,9),(9,10),(10,11),(11,12),(12,7), (13,14),(14,15),(15,16),(16,17),(17,18),(18,13), (19,20),(20,21),(21,22),(22,23),(23,24),(24,19)], []#=[(1,2,3,4,5,6), (7,8,9,10,11,12), (13,14,15,16,17,18), (19,20,21,22,23,24)]=#)
    p0 = [1 0 0; cos(2*pi/6) sin(2*pi/6) 0; cos(4*pi/6) sin(4*pi/6) 0; cos(6*pi/6) sin(6*pi/6) 0; cos(8*pi/6) sin(8*pi/6) 0; cos(10*pi/6) sin(10*pi/6) 0;
        1 0 0; cos(2*pi/6) 0 sin(2*pi/6); cos(4*pi/6) 0 sin(4*pi/6); cos(6*pi/6) 0 sin(6*pi/6); cos(8*pi/6) 0 sin(8*pi/6); cos(10*pi/6) 0 sin(10*pi/6);
        cos(10*pi/6) 0 sin(10*pi/6); 0 sqrt(1/2) -sqrt(1/2); cos(4*pi/6) sin(4*pi/6) 0; cos(4*pi/6) 0 sin(4*pi/6); 0 -sqrt(1/2) sqrt(1/2); cos(10*pi/6) sin(10*pi/6) 0;
        cos(8*pi/6) 0 sin(8*pi/6); 0 sqrt(1/2) -sqrt(1/2); cos(2*pi/6) sin(2*pi/6) 0; cos(2*pi/6) 0 sin(2*pi/6); 0 -sqrt(1/2) sqrt(1/2); cos(8*pi/6) sin(8*pi/6) 0]'
    initialConfiguration = toArray(p0,Weave)

    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave; colorscheme=[:yellow, :green, :cyan, :purple])
    toSphere(q, Weave)
end

function test_sphere_b2()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false, false,true,false,true,false,true, true,false,true,false,true,false, false,true,false,true,false,true], [(4, 10), (14, 20), (6, 24), (12, 19), (2, 15), (5, 18), (11, 13), (17, 23), (3, 21), (1, 7), (8, 16), (9, 22)], [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1), (7,8),(8,9),(9,10),(10,11),(11,12),(12,7), (13,14),(14,15),(15,16),(16,17),(17,18),(18,13), (19,20),(20,21),(21,22),(22,23),(23,24),(24,19)], []#=[(1,2,3,4,5,6), (7,8,9,10,11,12), (13,14,15,16,17,18), (19,20,21,22,23,24)]=#)
    p0 = [1 0 0; cos(2*pi/6) sin(2*pi/6) 0; cos(4*pi/6) sin(4*pi/6) 0; cos(6*pi/6) sin(6*pi/6) 0; cos(8*pi/6) sin(8*pi/6) 0; cos(10*pi/6) sin(10*pi/6) 0;
        1 0 0; cos(2*pi/6) 0 sin(2*pi/6); cos(4*pi/6) 0 sin(4*pi/6); cos(6*pi/6) 0 sin(6*pi/6); cos(8*pi/6) 0 sin(8*pi/6); cos(10*pi/6) 0 sin(10*pi/6);
        cos(8*pi/6) 0 sin(8*pi/6); 0 sqrt(1/2) -sqrt(1/2); cos(2*pi/6) sin(2*pi/6) 0; cos(2*pi/6) 0 sin(2*pi/6); 0 -sqrt(1/2) sqrt(1/2); cos(8*pi/6) sin(8*pi/6) 0;
        cos(10*pi/6) 0 sin(10*pi/6); 0 sqrt(1/2) -sqrt(1/2); cos(4*pi/6) sin(4*pi/6) 0; cos(4*pi/6) 0 sin(4*pi/6); 0 -sqrt(1/2) sqrt(1/2); cos(10*pi/6) sin(10*pi/6) 0;]'
    initialConfiguration = toArray(p0,Weave)
    #=bars = []  
    for i in 1:size(p0)[2], j in i+1:size(p0)[2]
        if norm(p0[:,i]-p0[:,j])<=10^(-10)
            push!(bars,(i,j))
        end
    end 
    println(collect(Set(bars)))=#

    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave; colorscheme=[:yellow, :green, :cyan, :purple])
end

function test_sphere_c()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true,false,true, false,true,false,true,false,true,false,true,false,true, true,false,true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true,false,true], [(5, 59), (14, 23), (19, 28), (25, 47), (34, 45), (41, 51), (2, 43), (35, 55), (24, 57), (3, 21), (13, 33), (12, 44), (6, 16), (39, 50), (27, 37), (29, 52), (1, 11), (18, 38), (7, 48), (20, 53), (9, 36), (17, 49), (40, 60), (30, 42), (10, 54), (22, 32), (15, 58), (46, 56), (8, 26), (4, 31)], 
                                            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(10,1), (11,12),(12,13),(13,14),(14,15),(15,16),(16,17),(17,18),(18,19),(19,20),(20,11), (21,22),(22,23),(23,24),(24,25),(25,26),(26,27),(27,28),(28,29),(29,30),(30,21), (31,32),(32,33),(33,34),(34,35),(35,36),(36,37),(37,38),(38,39),(39,40),(40,31), (41,42),(42,43),(43,44),(44,45),(45,46),(46,47),(47,48),(48,49),(49,50),(50,41), (51,52),(52,53),(53,54),(54,55),(55,56),(56,57),(57,58),(58,59),(59,60),(60,51)], 
                                            []#=[(1,2,3,4,5,6,7,8,9,10), (11,12,13,14,15,16,17,18,19,20), (21,22,23,24,25,26,27,28,29,30), (31,32,33,34,35,36,37,38,39,40), (41,42,43,44,45,46,47,48,49,50), (51,52,53,54,55,56,57,58,59,60)]=#; offset=0.1)
    p0 = [1 0 0; cos(2*pi/10) sin(2*pi/10) 0; cos(4*pi/10) sin(4*pi/10) 0; cos(6*pi/10) sin(6*pi/10) 0; cos(8*pi/10) sin(8*pi/10) 0; cos(10*pi/10) sin(10*pi/10) 0; cos(12*pi/10) sin(12*pi/10) 0; cos(14*pi/10) sin(14*pi/10) 0; cos(16*pi/10) sin(16*pi/10) 0; cos(18*pi/10) sin(18*pi/10) 0;
          1 0 0; cos(2*pi/10) 0 sin(2*pi/10); cos(4*pi/10) 0 sin(4*pi/10); cos(6*pi/10) 0 sin(6*pi/10); cos(8*pi/10) 0 sin(8*pi/10); cos(10*pi/10) 0 sin(10*pi/10); cos(12*pi/10) 0 sin(12*pi/10); cos(14*pi/10) 0 sin(14*pi/10); cos(16*pi/10) 0 sin(16*pi/10); cos(18*pi/10) 0 sin(18*pi/10);
          cos(4*pi/10) sin(4*pi/10) 0; 0 sqrt(1/2) sqrt(1/2); cos(6*pi/10) 0 sin(6*pi/10); -1/2 -sqrt(1/8) sqrt(5/8); -0.4 -sqrt(0.5) sqrt(0.34); cos(14*pi/10) sin(14*pi/10) 0; 0 -sqrt(1/2) -sqrt(1/2); cos(16*pi/10) 0 sin(16*pi/10); 1/2 sqrt(1/8) -sqrt(5/8); 1/2 sqrt(3/8) -sqrt(3/8);
          cos(6*pi/10) sin(6*pi/10) 0; 0 sqrt(1/2) sqrt(1/2); cos(4*pi/10) 0 sin(4*pi/10); 1/2 -sqrt(1/8) sqrt(5/8); 0.4 -sqrt(0.5) sqrt(0.34); cos(16*pi/10) sin(16*pi/10) 0; 0 -sqrt(1/2) -sqrt(1/2); cos(14*pi/10) 0 sin(14*pi/10); -1/2 sqrt(1/8) -sqrt(5/8); -1/2 sqrt(3/8) -sqrt(3/8);
          0 sqrt(3/8) -sqrt(5/8); 1/2 sqrt(3/8) -sqrt(3/8); cos(2*pi/10) sin(2*pi/10) 0; cos(2*pi/10) 0 sin(2*pi/10); 1/2 -sqrt(1/8) sqrt(5/8); 0 -sqrt(3/8) sqrt(5/8); -0.4 -sqrt(0.5) sqrt(0.34); cos(12*pi/10) sin(12*pi/10) 0; cos(12*pi/10) 0 sin(12*pi/10); -1/2 sqrt(1/8) -sqrt(5/8);
          0 sqrt(3/8) -sqrt(5/8); 1/2 sqrt(1/8) -sqrt(5/8); cos(18*pi/10) 0 sin(18*pi/10); cos(18*pi/10) sin(18*pi/10) 0; 0.4 -sqrt(0.5) sqrt(0.34); 0 -sqrt(3/8) sqrt(5/8); -1/2 -sqrt(1/8) sqrt(5/8); cos(8*pi/10) 0 sin(8*pi/10); cos(8*pi/10) sin(8*pi/10) 0; -1/2 sqrt(3/8) -sqrt(3/8)
    ]'

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:green,:teal,:red])
    toSphere(q, Weave)
end

function test_flattorus_plain()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true, false,true,false, true,false,true, false,true,false, true,false,true, false,true,false], [(1,10),(2,13),(3,16), (4,11),(5,14),(6,17), (7,12),(8,15),(9,18)], 
                                            [(1,2),(2,3), (4,5),(5,6), (7,8),(8,9), (10,11),(11,12), (13,14),(14,15), (16,17),(17,18)], 
                                            []; offset=0.1, manifold = "flattorus", samples = [[0.,0.,0],[]], PBC=[(3,1,[1.,0,0]), (6,4,[1.,0,0]), (9,7,[1.,0,0]), (12,10,[0.,1.,0]), (15,13,[0.,1.,0]), (18,16,[0.,1.,0])])
    p0 = [0 0 0; 0.5 0 0; 1 0 0; 0 0.5 0; 0.5 0.5 0; 1 0.5 0; 0 1 0; 0.5 1 0; 1 1 0; 
        0 0 0; 0 0.5 0; 0 1 0; 0.5 0 0; 0.5 0.5 0; 0.5 1 0; 1 0 0; 1 0.5 0; 1 1 0; 
    ]'

    initialConfiguration = toArray(p0,Weave)

    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave)

    q = computeOptimalWeaving(q, Weave; maxseconds=50)
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

function computeOptimalWeaving(initialConfiguration::Vector{Float64}, Weave::WeavingOnManifold; tol = 1e-3, maxseconds=250)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints)
    Q = x->energyFunction(x, Weave)
    G = ConstraintVariety(Weave.coordinateVariables, Weave.constraints, length(Weave.coordinateVariables), length(Weave.coordinateVariables)-length(Weave.constraints))
    @time resultmin = findminima(q, 1e-3, G, Q; whichstep="gaussnewtonstep", stepdirection = "gradientdescent", maxseconds=maxseconds)
    @time resultmin = findminima(q, 1e-3, G, Q; whichstep="EDStep", homotopyMethod="HomotopyContinuation", stepdirection = "gradientdescent", maxseconds=maxseconds)
    return(resultmin.computedpoints[end])
end

test_sphere_a()

#TODO how does the energy behave, when I squish the sphere => become an ellipse. 
#TODO how to restrict weavings on manifolds?
#TODO Tensegrities in the plane with saddle?
#TODO What if the cables are under compression instead of tension?W
#TODO What is the optimal manifold for the weaving => random 4-coordinated graph, alternating. Fixed volume.
#TODO Hyperauxeticity paper
end
