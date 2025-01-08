module WeavingsOnManifolds
import LinearAlgebra: pinv, norm, cross, nullspace, I
#import GLMakie: hidedecorations!, hidespines!, Figure, Axis3, hidespines!, hidedecorations!, Point3f0, scatter!, linesegments!, RGBA, save, Axis, lines!, axislegend
#import Implicit3DPlotting: plot_implicit_surface!
import Plots

include("../../HomotopyOpt.jl/src/HomotopyOpt.jl")
import LaTeXStrings: @L_str
import HomotopyContinuation: Expression, Variable, @var, evaluate, differentiate, System, real_solutions, solve
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
    cablejoints
    pinned_vertices
    ellipse
    struts

    function WeavingOnManifold(offsetList::Vector, bars::Vector, cables::Vector, planes::Vector; rotationsmatrix=Matrix{Float64}(I,3,3), manifold="sphere", ellipse=[1,1], implicitManifold=x->x[1]^2+x[2]^2+x[3]^2-1, offset=0.1, torusknot=false, angleConstraints=false, samples=[], PBC=[], pinned_vertices=Dict())
        (offset<=0.5 && offset>=0.05) || @warn "We recommend an offset between 5% and 50%."
        manifold!="free" || typeof(implicitManifold([1.,2,3]))==Float64
        ([i for i in 1:length(offsetList)] == sort(vcat([bar[1] for bar in bars], [bar[2] for bar in bars])) && (Set(1:length(offsetList)) == Set(vcat([cable[1] for cable in cables], [cable[2] for cable in cables]))))# || throw(error("Cables and Bars don't cover every vertex"))

        @var x[1:3, 1:length(offsetList)]

        xs = zeros(Expression,(3,length(offsetList)))
        manifoldEquations, barEquations, angleEquations, planeEquations, periodicBoundaryConditions, xvarz, positions, outsideindices, cableJoints, struts = Vector{Expression}([]), Vector{Expression}([]), Vector{Expression}([]), Vector{Expression}([]), Vector{Expression}([]), Vector{Variable}([]), [], [], [], []

        for i in findall(t->t==1,offsetList), j in findall(t->t==1,offsetList)
            if j<=i
                continue
            end
            push!(outsideindices, (i,j))
        end

        for cable in cables
            strut = (offsetList[cable[1]] ? cable[1] : (bars[findfirst(bar->cable[1] in bar, bars)][1]==cable[1] ? bars[findfirst(bar->cable[1] in bar, bars)][2] : bars[findfirst(bar->cable[1] in bar, bars)][1]), offsetList[cable[2]] ? cable[2] : (bars[findfirst(bar->cable[2] in bar, bars)][1]==cable[2] ? bars[findfirst(bar->cable[2] in bar, bars)][2] : bars[findfirst(bar->cable[2] in bar, bars)][1]))
            push!(struts, strut[1]<strut[2] ? strut : (strut[2],strut[1]))
        end
        struts = collect(Set(struts))
        display(struts)

        if isequal(lowercase(manifold), "sphere")
            #xs[:,1] = [offsetList[1] ? 1+offset : 1-offset, 0, 0]
            #xs[:,2] = [0, x[2,2], x[3,2]]

            xs[:,1] = isempty(samples) ? (1+offset*(2*offsetList[1]-1)) * [ 1,0,0] : (1+offset*(2*offsetList[1]-1)) .* samples[1] ./ norm(samples[1])
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
        elseif isequal(lowercase(manifold), "ellipsoid")
            q0 = isempty(samples) ? (1+offset*(2*offsetList[1]-1)) * rotationsmatrix * [sqrt(ellipse[1]),0,0] : samples[1]
            println(q0)
            @var start_var[1:3]
            q0 = newtonCorrect(q0, start_var, [sum([1/sqrt(ellipse[1]),1/sqrt(ellipse[1]),ellipse[1]] .* (start_var).^2) - (1+offset*(2*offsetList[1]-1))^2])
            xs[1,1] = q0[1]
            xs[2:3,1] = x[2:3,1]
            xs[:,2:length(offsetList)] .= x[:,2:length(offsetList)]
    
            for bar in bars
                xs[:,bar[2]] = offsetList[bar[2]] ? xs[:,bar[1]] + 2*offset/sqrt((1/(sqrt(ellipse[1]))*xs[1,bar[1]])^2+(1/(sqrt(ellipse[1]))*xs[2,bar[1]])^2+(ellipse[1]*xs[3,bar[1]])^2) * [1/sqrt(ellipse[1]),1/sqrt(ellipse[1]),ellipse[1]] .* xs[:,bar[1]] : xs[:,bar[1]] - 2*offset/sqrt((1/(sqrt(ellipse[1]))*xs[1,bar[1]])^2+(1/(sqrt(ellipse[1]))*xs[2,bar[1]])^2+(ellipse[1]*xs[3,bar[1]])^2) * [1/sqrt(ellipse[1]),1/sqrt(ellipse[1]),ellipse[1]] .* xs[:,bar[1]]
            end
            for i in 1:size(xs)[2], j in 1:size(xs)[1]
                try
                    push!(xvarz, Variable(xs[j,i]))
                    push!(positions, (j,i))
                catch
                    continue
                end
            end
    
            samples = [q0]
            manifoldEquations = [offsetList[bar[1]] ? sum([1/sqrt(ellipse[1]),1/sqrt(ellipse[1]),ellipse[1]] .* (xs[:,bar[1]]).^2)-(1+offset)^2 : sum([1/sqrt(ellipse[1]),1/sqrt(ellipse[1]),ellipse[1]] .* (xs[:,bar[1]]).^2)-(1-offset)^2 for bar in bars]
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
            xs[1:1,1] = isempty(samples) ? -(1.5+(2*offsetList[1]-1)*Weave.offset) : samples[1][1]
            #xs[:,2] = [torusknot[1] ? cos(2*pi/torusknot[2])*x[2,2] : 0, torusknot[1] ? sin(2*pi/torusknot[2])*x[2,2] : x[2,2], x[3,2]]
            xs[:,2:length(offsetList)] .= x[:,2:length(offsetList)]
            
            if !torusknot
                append!(barEquations, [sum((xs[:,bar[1]]-xs[:,bar[2]]).^2)-(2*offset)^2 for bar in bars])
            else
                append!(barEquations, [sum((xs[:,bar[1]]-xs[:,bar[2]]).^2)-(1)^2 for bar in bars])
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
        elseif isequal(lowercase(manifold), "cylinder")
            @var h
            xs[:,1] = isempty(samples) ? (1+offset*(2*offsetList[1]-1)) * [ 1,0,0] : vcat((1+offset*(2*offsetList[1]-1)) .* samples[1][1:2], samples[1][3])
            xs[:,2:length(offsetList)] .= x[:,2:length(offsetList)]
            for bar in bars
                if bar[1] in [entry[2] for entry in PBC]
                    cond = filter(entry->entry[2]==bar[1], PBC)[1]
                    xs[:,cond[2]] = xs[:,cond[1]] + h*cond[3]
                end
                xs[:,bar[2]] = offsetList[bar[2]]==1 ? [xs[1,bar[1]]*(1+offset)/(1-offset), xs[2,bar[1]]*(1+offset)/(1-offset), xs[3,bar[1]]] : [xs[1,bar[1]]*(1-offset)/(1+offset),xs[2,bar[1]]*(1-offset)/(1+offset),xs[3,bar[1]]]
            end
            for cond in PBC
                xs[:,cond[2]] = xs[:,cond[1]] + h*cond[3]
            end
            for i in 1:size(xs)[2], j in 1:size(xs)[1]
                try
                    #=if any(t->t==Variable(xs[j,i]), xvarz)
                        continue
                    end=#
                    push!(xvarz, Variable(xs[j,i]))
                    push!(positions, (j,i))
                catch
                    continue
                end
            end
            push!(xvarz,h)

            samples = isempty(samples) ? append!(samples, [(1+2*offset*(offsetList[1]-0.5)) * [ 1,0,0]]) : samples
            manifoldEquations = [sum( xs[1:2,bar[1]].^2 ) - (1+2*offset*(offsetList[bar[1]]-0.5))^2 for bar in bars]
            manifoldEquations = filter(t->t!=0, manifoldEquations)
            for key in keys(pinned_vertices)
                for i in 1:length(pinned_vertices[key][1])
                    push!(manifoldEquations, xs[Int(pinned_vertices[key][1][i]),key] - pinned_vertices[key][2][i])
                end
            end
        elseif isequal(lowercase(manifold), "hyperboloid")    
                #xs[1:3,1] = isempty(samples) ? (1+offset*(2*offsetList[1]-1)) * [ 1,0,0] : samples[1]
                xs[:,1:length(offsetList)] .= x[:,1:length(offsetList)]
    
                for key in keys(pinned_vertices)
                    xs[:,key] = pinned_vertices[key]
                end
                append!(barEquations, [sum((xs[:,bar[1]]-xs[:,bar[2]]).^2)-(2*offset)^2 for bar in bars])

                for i in 1:size(xs)[2], j in 1:size(xs)[1]
                    try
                        push!(xvarz, Variable(xs[j,i]))
                        push!(positions, (j,i))
                    catch
                        continue
                    end
                end
    
                samples = isempty(samples) ? append!(samples, [offsetList[1] ? (1+offset) * [ 1,0,0] : (1-offset) * [ 1,0,0]]) : samples
                manifoldEquations = [sum( xs[1:2,i].^2 ) - xs[3,i]^2 + ((1+offset*(2*offsetList[i]-1)))^2 for i in 2:size(xs)[2] if !(i in collect(keys(pinned_vertices)))]
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
        new(lowercase(manifold), offsetList, offset, filter(t->t!=0,Vector{Expression}(vcat(manifoldEquations, barEquations, planeEquations, angleEquations, periodicBoundaryConditions))), Vector{Variable}(xvarz), bars, cables, planes, positions, samples, outsideindices, PBC, cableJoints, pinned_vertices, ellipse, struts)
    end
end

function toMatrix(configuration, Weave::WeavingOnManifold)
    p0 = zeros(Number, 3, Int(length(Weave.offsetList)))
    global count = 1
    if Weave.manifold=="sphere"
        p0[:,1] = (1+(2*Weave.offsetList[1]-1)*Weave.offset) * Weave.samples[1] ./ norm(Weave.samples[1])
    elseif Weave.manifold=="cylinder"
        p0[:,1] = Weave.samples[1]
    elseif Weave.manifold=="torus"
        p0[1,1] = isempty(Weave.samples) ? -(1.5+(2*Weave.offsetList[1]-1)*Weave.offset) : Weave.samples[1][1]
    else
        p0[:,1] = Weave.samples[1]
    end
    for pos in Weave.variablePositions
        p0[pos[1],pos[2]] = configuration[count]
        global count = count+1
    end


    if Weave.manifold=="sphere"
        for bar in Weave.bars
            p0[:,bar[2]] = Weave.offsetList[bar[2]] ? p0[:,bar[1]]*(1+Weave.offset)/(1-Weave.offset) : p0[:,bar[1]]*(1-Weave.offset)/(1+Weave.offset)
        end
    elseif Weave.manifold=="ellipsoid"
        for bar in Weave.bars
            p0[:,bar[2]] = Weave.offsetList[bar[2]] ? p0[:,bar[1]] + 2*Weave.offset/sqrt(2*(1/sqrt(Weave.ellipse[1]))^2+Weave.ellipse[1]^2)*[1/sqrt(Weave.ellipse[1]),1/sqrt(Weave.ellipse[1]),Weave.ellipse[1]] .* p0[:,bar[1]] / norm(p0[:,bar[1]]) : p0[:,bar[1]] - 2*Weave.offset/sqrt((1/sqrt(Weave.ellipse[1]))^2+1/sqrt(Weave.ellipse[1])^2+Weave.ellipse[1]^2) * [1/sqrt(Weave.ellipse[1]),1/sqrt(Weave.ellipse[1]),Weave.ellipse[1]] .* p0[:,bar[1]] / norm(p0[:,bar[1]])
        end
    elseif Weave.manifold=="cylinder"
        for bar in Weave.bars
            if bar[1] in [entry[2] for entry in Weave.periodicBoundary]
                cond = filter(entry->entry[2]==bar[1], Weave.periodicBoundary)[1]
                p0[:,cond[2]] = p0[:,cond[1]] + configuration[end]*cond[3]
            end
            p0[:,bar[2]] = vcat((1+(2*Weave.offsetList[bar[2]]-1)*Weave.offset)/(1-(2*Weave.offsetList[bar[2]]-1)*Weave.offset) * p0[1:2,bar[1]], p0[3,bar[1]])
        end
        for cond in Weave.periodicBoundary
            p0[:,cond[2]] = p0[:,cond[1]]+configuration[end]*cond[3]
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
    if Weave.manifold=="cylinder"
        push!(configuration, norm(p[:,Weave.periodicBoundary[1][2]]-p[:,Weave.periodicBoundary[1][1]]))
    end
    return configuration
end

function energyFunction(configuration, Weave::WeavingOnManifold; no_struts=false, no_cables = false, bending = false)
    p = toMatrix(configuration, Weave)
    Q = 0
    if !no_cables
        Q += 1/2*sum([ sum( (p[:,cable[1]] - p[:,cable[2]]).^2 ) for cable in Weave.cables])
    end
    #Q += -2*sum([ sum( (p[:,t[1]] - p[:,t[2]]).^2 ) for t in Weave.outsideindices])
    if !no_struts
        Q += 8/length(configuration)*sum([ 1 / norm((p[:,t[1]] - p[:,t[2]])) for t in Weave.outsideindices])
        #Q += 0.5*sum([ 1 / norm((p[:,t[1]] - p[:,t[2]])) for t in Weave.struts])
    end
    if bending
        for edge in Weave.cables
            cable1 = filter(cable->edge[1] in cable && cable!=edge, Weave.cables)
            if !isempty(cable1)
                cable1 = cable1[1][2]==edge[1] ? cable1[1] : (cable1[1][2],cable1[1][1])
                Q += 0.5*(1-(p[:,cable1[2]] - p[:,cable1[1]])'*(p[:,edge[2]] - p[:,edge[1]])/(norm(p[:,cable1[2]] - p[:,cable1[1]])*norm(p[:,edge[2]] - p[:,edge[1]])))
            end
            cable2 = filter(cable->edge[2] in cable && cable!=edge, Weave.cables)
            if !isempty(cable2)
                cable2 = cable2[1][1]==edge[2] ? cable2[1] : (cable2[1][2],cable2[1][1])
                Q += 0.5*(1-(p[:,cable2[2]] - p[:,cable2[1]])'*(p[:,edge[2]] - p[:,edge[1]])/(norm(p[:,cable2[2]] - p[:,cable2[1]])*norm(p[:,edge[2]] - p[:,edge[1]])))
            end
        end
    end
    #TODO Thompson potential on the sphere vs. Lennard-Jones?
    return Q
end

function plotWeaving(configuration::Vector{Float64}, Weave::WeavingOnManifold; colorscheme=[])
    p0 = toMatrix(configuration, Weave)
    fig = Figure(size = (1000,1000))
    ax = Axis3(fig[1,1], aspect=(1.,1,1))
    hidedecorations!(ax)
    hidespines!(ax)
   # hidespines!(ax); hidedecorations!(ax);
    implicit_manifold = x->(Weave.manifold == "torus") ? (0.75 + (x[1]^2+x[2]^2+x[3]^2))^2 - 4*(x[1]^2+x[2]^2) : ((Weave.manifold == "flattorus") ? x[3] : x[1]^2+x[2]^2+x[3]^2-0.95)
    plot_implicit_surface!(ax, implicit_manifold; wireframe=false, transparency=true, color=RGBA(0.85,0.85,0.85,0.6), samples = (75,75,75), xlims = (Weave.manifold == "flattorus") ? (0,1) : (-1-Weave.offset, 1+Weave.offset), ylims = (Weave.manifold == "flattorus") ? (0,1) : (-1-Weave.offset, 1+Weave.offset), zlims = (Weave.manifold == "flattorus") ? (-2*Weave.offset-0.01,2*Weave.offset+0.01) : (-1-Weave.offset, 1+Weave.offset))
    foreach(bar->linesegments!(ax, [Point3f0(p0[:,bar[1]]), Point3f0(p0[:,bar[2]])]; color=:grey40, linewidth=13), Weave.bars)

    if isempty(colorscheme)
        foreach(cable->linesegments!(ax, [Point3f0(p0[:,cable[1]]), Point3f0(p0[:,cable[2]])]; color=:red, linewidth=13), Weave.cables)
    else
        cableJoints = []
        for cable in Weave.cables
            if !any(t->(cable[1] in t)||(cable[2] in t), cableJoints)
                push!(cableJoints, [cable[1],cable[2]])
            else
                append!(cableJoints[findfirst(t->(cable[1] in t)||(cable[2] in t), cableJoints)], [cable[1],cable[2]])
            end
        end
        foreach(cable->linesegments!(ax, [Point3f0(p0[:,cable[1]]), Point3f0(p0[:,cable[2]])]; color=colorscheme[findfirst(t->cable[1] in t, cableJoints)], linewidth=13), Weave.cables)
    end

    scatter!(ax, [Point3f0(p0[:,i]) for i in 1:size(p0)[2]]; color=:black, markersize=35)
    #save("currentweaving.png",fig)
    display(fig)
end

function doubleRings(configuration, Weave)
    offsetList_ = []#vcat(Weave.offsetList, [!offset for offset in Weave.offsetList], [!offset for offset in Weave.offsetList], Weave.offsetList)
    bars_ = []#vcat(Weave.bars, [(bar[1]+length(Weave.offsetList), bar[2]+length(Weave.offsetList)) for bar in Weave.bars], [(bar[1]+2*length(Weave.offsetList), bar[2]+2*length(Weave.offsetList)) for bar in Weave.bars], [(bar[1]+3*length(Weave.offsetList), bar[2]+3*length(Weave.offsetList)) for bar in Weave.bars])
    cables_ = []#vcat(Weave.cables, [(cable[1]+length(Weave.offsetList), cable[2]+length(Weave.cables)) for cable in Weave.cables], [(cable[1]+2*length(Weave.offsetList), bar[2]+2*length(Weave.offsetList)) for bar in Weave.bars], [(bar[1]+3*length(Weave.offsetList), bar[2]+3*length(Weave.offsetList)) for bar in Weave.bars])
    configuration_ = []
    for bar in Weave.bars
        adjacent_cables = vcat(findall(cable->bar[1] in cable, Weave.cables), findall(cable->bar[2] in cable, Weave.cables))
        normals = [cross(configuration[:,Weave.cables[adjacent_cables[1]][1]]-configuration[:,Weave.cables[adjacent_cables[1]][2]], configuration[:,Weave.cables[adjacent_cables[2]][1]]-configuration[:,Weave.cables[adjacent_cables[2]][2]]),
                    cross(configuration[:,Weave.cables[adjacent_cables[3]][1]]-configuration[:,Weave.cables[adjacent_cables[3]][2]], configuration[:,Weave.cables[adjacent_cables[4]][1]]-configuration[:,Weave.cables[adjacent_cables[4]][2]])]
        normals = [normal ./norm(normal) for normal in normals]
        
        append!(configuration_, vcat([configuration[:,bar[1]]+0.3*(normals[1]+normals[2]), configuration[:,bar[1]]+0.3*(normals[1]-normals[2]), configuration[:,bar[1]]+0.3*(-normals[1]+normals[2]), configuration[:,bar[1]]+0.3*(-normals[1]-normals[2])],
                                        [configuration[:,bar[2]]+0.3*(normals[1]+normals[2]), configuration[:,bar[2]]+0.3*(normals[1]-normals[2]), configuration[:,bar[2]]+0.3*(-normals[1]+normals[2]), configuration[:,bar[2]]+0.3*(-normals[1]-normals[2])]))
        append!(bars_,[(i,i+4) for i in length(configuration_)-7:length(configuration_)-4])
        append!(offsetList_,[Weave.offsetList[bar[1]], Weave.offsetList[bar[1]], Weave.offsetList[bar[1]], Weave.offsetList[bar[1]],
                        !Weave.offsetList[bar[2]], !Weave.offsetList[bar[2]], !Weave.offsetList[bar[2]], !Weave.offsetList[bar[2]]])
        if abs(normals[1]'*(configuration_[length(configuration_)-7]-configuration_[length(configuration_)-2])) < abs(normals[2]'*(configuration_[length(configuration_)-7]-configuration_[length(configuration_)-2]))
            append!(cables_, [(length(configuration_)-7,length(configuration_)-2), (length(configuration_)-4, length(configuration_)-1), (length(configuration_)-5,length(configuration_)-3), (length(configuration_)-6,length(configuration_))])
        else
            append!(cables_, [(length(configuration_)-7,length(configuration_)-1), (length(configuration_)-4, length(configuration_)-2), (length(configuration_)-5,length(configuration_)), (length(configuration_)-6,length(configuration_)-3)])
        end
    end
    Weave_ = WeavingsOnManifolds.WeavingOnManifold(offsetList_, bars_, cables_, [])
    configuration_= Matrix{Float64}(hcat(configuration_...))
    println(offsetList_)
    println(cables_)
    println(bars_)
    println(configuration_)
    return configuration_
end

function toSphere(configuration, Weave; pltstring="newconf")
    p0 = toMatrix(configuration, Weave)
    newEdgeList = []
    newPointList = []
    outer = 1
    #epsilon = 0.15
    for edge in Weave.cables
        outer = norm(p0[:,edge[1]])>1 ? 1 : -1
        cable1, cable2 = filter(cable->edge[1] in cable && cable!=edge, Weave.cables)[1], filter(cable->edge[2] in cable && cable!=edge, Weave.cables)[1]
        v1, v2 = p0[:,edge[2]]-p0[:,edge[1]] + p0[:,cable1[2]]-p0[:,cable1[1]], p0[:,cable2[2]]-p0[:,cable2[1]] + p0[:,edge[2]]-p0[:,edge[1]] 
        v1, v2 = v1 - p0[:,edge[1]] * (p0[:,edge[1]]' * v1), v2 - p0[:,edge[2]] * (p0[:,edge[2]]' * v2)
        v1 = 0.65*v1; v2=0.65*v2
        A = [([3 2 1 0; 0 0 1 0; 1 1 1 1; 0 0 0 1] \ [v2[i],v1[i],p0[:,edge[2]][i],p0[:,edge[1]][i]]) for i in 1:3]

        for t in 0:0.05:0.95
            curP = t*(p0[:,edge[2]]-p0[:,edge[1]])+p0[:,edge[1]]
            #push!(newPointList, (outer==1 ? (1+epsilon-(6*epsilon)*t^2+(4*epsilon)*t^3) : (1+epsilon-(6*epsilon)*(1-t)^2+(4*epsilon)*(1-t)^3))*curP./norm(curP))
            push!(newPointList, [A[i]'*[t^3,t^2,t,1] for i in 1:length(A)])
            #push!(newPointList, curP ./ norm(curP))
            push!(newEdgeList, (length(newPointList), length(newPointList)+1))
        end
        curP = p0[:,edge[2]]
        push!(newPointList, p0[:,edge[2]])#(1-outer*epsilon)*curP)
    end
    for i in size(p0)[2]
        push!(newPointList, p0[:,i])# ./ norm(p0[:,i]))
    end
    open("weavingframework$(pltstring).poly", "w") do f
        write(f, "POINTS\n")
        foreach(i->write(f, string("$(i): ", newPointList[i][1], " ", newPointList[i][2], " ", newPointList[i][3],"\n")), 1:length(newPointList))
        write(f,"POLYS\n")
        foreach(i->write(f, string("$(i): ", newEdgeList[i][1], " ", newEdgeList[i][2],"\n")), 1:length(newEdgeList))
        write(f,"END")
    end
end


function toEllipsoid(configuration, Weave; ellipse, pltstring = "ellispoid_oct1.0")
    p0 = toMatrix(configuration, Weave)
    newEdgeList = []
    newPointList = []
    for edge in Weave.cables
        cable1, cable2 = filter(cable->edge[1] in cable && cable!=edge, Weave.cables)[1], filter(cable->edge[2] in cable && cable!=edge, Weave.cables)[1]
        v1, v2 = p0[:,edge[2]]-p0[:,edge[1]] + p0[:,cable1[2]]-p0[:,cable1[1]], p0[:,cable2[2]]-p0[:,cable2[1]] + p0[:,edge[2]]-p0[:,edge[1]] 
        v1, v2 = v1 - 1/((1/(ellipse[1]*ellipse[2]))^2+ellipse[1]^2+ellipse[2]^2)*([1/(ellipse[1]*ellipse[2]),ellipse[1],ellipse[2]] .* p0[:,edge[1]]) ./ norm(p0[:,edge[1]]) * (([1/(ellipse[1]*ellipse[2]),ellipse[1],ellipse[2]] .* p0[:,edge[1]] ./ norm(p0[:,edge[1]]))' * v1), v2 - 1/((1/(ellipse[1]*ellipse[2]))^2+ellipse[1]^2+ellipse[2]^2)*([1/(ellipse[1]*ellipse[2]),ellipse[1],ellipse[2]] .* p0[:,edge[2]] ./ norm(p0[:,edge[1]])) * (([1/(ellipse[1]*ellipse[2]),ellipse[1],ellipse[2]] .* p0[:,edge[2]] ./ norm(p0[:,edge[1]]))' * v2)
        v1 = 0.7*v1; v2=0.7*v2
        A = [([3 2 1 0; 0 0 1 0; 1 1 1 1; 0 0 0 1] \ [v2[i],v1[i],p0[:,edge[2]][i],p0[:,edge[1]][i]]) for i in 1:3]

        for t in 0:0.05:0.95
            curP = t*(p0[:,edge[2]]-p0[:,edge[1]])+p0[:,edge[1]]
            #push!(newPointList, (outer==1 ? (1+epsilon-(6*epsilon)*t^2+(4*epsilon)*t^3) : (1+epsilon-(6*epsilon)*(1-t)^2+(4*epsilon)*(1-t)^3))*curP./norm(curP))
            push!(newPointList, [A[i]'*[t^3,t^2,t,1] for i in 1:length(A)])
            #push!(newPointList, curP ./ norm(curP))
            push!(newEdgeList, (length(newPointList), length(newPointList)+1))
        end
        curP = p0[:,edge[2]]
        push!(newPointList, p0[:,edge[2]])#(1-outer*epsilon)*curP)
    end
    for i in size(p0)[2]
        push!(newPointList, p0[:,i])# ./ norm(p0[:,i]))
    end
    open("weavingframework$(pltstring).poly", "w") do f
        write(f, "POINTS\n")
        foreach(i->write(f, string("$(i): ", newPointList[i][1], " ", newPointList[i][2], " ", newPointList[i][3],"\n")), 1:length(newPointList))
        write(f,"POLYS\n")
        foreach(i->write(f, string("$(i): ", newEdgeList[i][1], " ", newEdgeList[i][2],"\n")), 1:length(newEdgeList))
        write(f,"END")
    end
end


function toHyperboloid(configuration, Weave)
    p0 = toMatrix(configuration, Weave)
    newEdgeList = []
    newPointList = []
    outer = 1
    @var y[1:3]
    hyperboloid_normal = differentiate(y[1]^2+y[2]^2-y[3]^2+1, y)
    epsilon = 0.15
    for edge in Weave.cables
        outer = norm(p0[:,edge[1]])>1 ? 1 : -1
        cable1, cable2 = filter(cable->edge[1] in cable && cable!=edge, Weave.cables), filter(cable->edge[2] in cable && cable!=edge, Weave.cables)
        if isempty(cable1)
            v1 = p0[:,edge[2]]-p0[:,edge[1]]
            n1 = [p0[1,edge[1]], p0[2,edge[1]], 0]
        else
            v1 = p0[:,edge[2]]-p0[:,edge[1]] + p0[:,cable1[1][2]]-p0[:,cable1[1][1]]
            n1 = evaluate(hyperboloid_normal, y=>p0[:,edge[1]])
        end

        if isempty(cable2)
            v2 = p0[:,edge[2]]-p0[:,edge[1]]
            n2 = [p0[1,edge[2]], p0[2,edge[2]], 0]
        else
            v2 = p0[:,cable2[1][2]]-p0[:,cable2[1][1]] + p0[:,edge[2]]-p0[:,edge[1]] 
            n2 = evaluate(hyperboloid_normal, y=>p0[:,edge[2]])
        end

        n1, n2 = n1 ./ norm(n1), n2 ./ norm(n2)
        v1, v2 = v1 - n1 * (n1' * v1), v2 - n2 * (n2' * v2)
        v1 = 0.5*v1; v2=0.5*v2
        A = [([3 2 1 0; 0 0 1 0; 1 1 1 1; 0 0 0 1] \ [v2[i],v1[i],p0[:,edge[2]][i],p0[:,edge[1]][i]]) for i in 1:3]

        for t in 0:0.05:0.95
            push!(newPointList, [A[i]'*[t^3,t^2,t,1] for i in 1:length(A)])
            push!(newEdgeList, (length(newPointList), length(newPointList)+1))
        end
        curP = p0[:,edge[2]]
        push!(newPointList, p0[:,edge[2]])#(1-outer*epsilon)*curP)
    end
    for i in size(p0)[2]
        push!(newPointList, p0[:,i])# ./ norm(p0[:,i]))
    end
    open("weavingframeworkmodelonhyperboloid.poly", "w") do f
        write(f, "POINTS\n")
        foreach(i->write(f, string("$(i): ", newPointList[i][1], " ", newPointList[i][2], " ", newPointList[i][3],"\n")), 1:length(newPointList))
        write(f,"POLYS\n")
        foreach(i->write(f, string("$(i): ", newEdgeList[i][1], " ", newEdgeList[i][2],"\n")), 1:length(newEdgeList))
        write(f,"END")
    end
end

function toTorus(configuration, Weave)
    p0 = toMatrix(configuration, Weave)
    newEdgeList = []
    newPointList = []
    for i in 1:size(p0)[2]
        #p0[:,i] = p0[:,i]./norm(p0[:,i])
    end
    outer = 1
    epsilon = 0.15
    @var y[1:3]
    torus_normal = differentiate(((1-(0.5)^2)+sum(y.^2))^2-4*sum(y[1:2].^2), y)
    
    for edge in Weave.cables
        outer = norm(p0[:,edge[1]])>1 ? 1 : -1
        cable1, cable2 = filter(cable->edge[1] in cable && cable!=edge, Weave.cables)[1], filter(cable->edge[2] in cable && cable!=edge, Weave.cables)[1]
        v1, v2 = p0[:,edge[2]]-p0[:,edge[1]] + p0[:,cable1[2]]-p0[:,cable1[1]], p0[:,cable2[2]]-p0[:,cable2[1]] + p0[:,edge[2]]-p0[:,edge[1]] 
        n1, n2 = evaluate(torus_normal, y=>p0[:,edge[1]]), evaluate(torus_normal, y=>p0[:,edge[2]])
        n1, n2 = n1 ./ norm(n1), n2 ./ norm(n2)
        v1, v2 = v1 - n1 * (n1' * v1), v2 - n2 * (n2' * v2)
        v1 = .6*v1; v2= 0.6*v2
        A = [([3 2 1 0; 0 0 1 0; 1 1 1 1; 0 0 0 1] \ [v2[i],v1[i],p0[:,edge[2]][i],p0[:,edge[1]][i]]) for i in 1:3]

        for t in 0:0.05:0.95
            curP = t*(p0[:,edge[2]]-p0[:,edge[1]])+p0[:,edge[1]]
            #push!(newPointList, (outer==1 ? (1+epsilon-(6*epsilon)*t^2+(4*epsilon)*t^3) : (1+epsilon-(6*epsilon)*(1-t)^2+(4*epsilon)*(1-t)^3))*curP./norm(curP))
            push!(newPointList, [A[i]'*[t^3,t^2,t,1] for i in 1:length(A)])
            #push!(newPointList, curP ./ norm(curP))
            push!(newEdgeList, (length(newPointList), length(newPointList)+1))
        end
        curP = p0[:,edge[2]]
        push!(newPointList, p0[:,edge[2]])#(1-outer*epsilon)*curP)
    end
    for i in size(p0)[2]
        push!(newPointList, p0[:,i])# ./ norm(p0[:,i]))
    end
    open("weavingframeworkmodelontorustrefoil.poly", "w") do f
        write(f, "POINTS\n")
        foreach(i->write(f, string("$(i): ", newPointList[i][1], " ", newPointList[i][2], " ", newPointList[i][3],"\n")), 1:length(newPointList))
        write(f,"POLYS\n")
        foreach(i->write(f, string("$(i): ", newEdgeList[i][1], " ", newEdgeList[i][2],"\n")), 1:length(newEdgeList))
        write(f,"END")
    end
end

function toFlatTorus(configuration, Weave)
    p0 = toMatrix(configuration, Weave)
    newEdgeList = []
    newPointList = []
    for i in 1:size(p0)[2]
        #p0[:,i] = p0[:,i]./norm(p0[:,i])
    end

    outer = 1
    epsilon = 0.15
    @var y[1:3]
    torus_normal = differentiate(((1-(0.5)^2)+sum(y.^2))^2-4*sum(y[1:2].^2), y)
    
    for edge in Weave.cables
        outer = norm(p0[:,edge[1]])>1 ? 1 : -1
        cable1, cable2 = filter(cable->edge[1] in cable && cable!=edge, Weave.cables), filter(cable->edge[2] in cable && cable!=edge, Weave.cables)
        cable1, cable2 = (isempty(cable1) ? cable2[1] : cable1[1]), (isempty(cable2) ? cable1[1] : cable2[1])
        v1, v2 = p0[:,edge[2]]-p0[:,edge[1]] + p0[:,cable1[2]]-p0[:,cable1[1]], p0[:,cable2[2]]-p0[:,cable2[1]] + p0[:,edge[2]]-p0[:,edge[1]] 
        n1, n2 = [0,0,1], [0,0,1]
        v1, v2 = v1 - n1 * (n1' * v1), v2 - n2 * (n2' * v2)
        v1 = .6*v1; v2= 0.6*v2
        A = [([3 2 1 0; 0 0 1 0; 1 1 1 1; 0 0 0 1] \ [v2[i],v1[i],p0[:,edge[2]][i],p0[:,edge[1]][i]]) for i in 1:3]

        for t in 0:0.05:0.95
            curP = t*(p0[:,edge[2]]-p0[:,edge[1]])+p0[:,edge[1]]
            #push!(newPointList, (outer==1 ? (1+epsilon-(6*epsilon)*t^2+(4*epsilon)*t^3) : (1+epsilon-(6*epsilon)*(1-t)^2+(4*epsilon)*(1-t)^3))*curP./norm(curP))
            push!(newPointList, [A[i]'*[t^3,t^2,t,1] for i in 1:length(A)])
            #push!(newPointList, curP ./ norm(curP))
            push!(newEdgeList, (length(newPointList), length(newPointList)+1))
        end
        curP = p0[:,edge[2]]
        push!(newPointList, p0[:,edge[2]])#(1-outer*epsilon)*curP)
    end
    for i in size(p0)[2]
        push!(newPointList, p0[:,i])# ./ norm(p0[:,i]))
    end
    open("weavingframeworkmodelontorustrefoil.poly", "w") do f
        write(f, "POINTS\n")
        foreach(i->write(f, string("$(i): ", newPointList[i][1], " ", newPointList[i][2], " ", newPointList[i][3],"\n")), 1:length(newPointList))
        write(f,"POLYS\n")
        foreach(i->write(f, string("$(i): ", newEdgeList[i][1], " ", newEdgeList[i][2],"\n")), 1:length(newEdgeList))
        write(f,"END")
    end
end

function toCylinder(configuration, Weave)
    p0 = toMatrix(configuration, Weave)
    newEdgeList = []
    newPointList = []
    outer = 1
    epsilon = 0.125
    @var y[1:3]
    
    for edge in Weave.cables
        outer = norm(p0[:,edge[1]])>1 ? 1 : -1
        cable1, cable2 = filter(cable->edge[1] in cable && cable!=edge, Weave.cables), filter(cable->edge[2] in cable && cable!=edge, Weave.cables)
        if isempty(cable1) 
            if filter(t->edge[1] in t, Weave.periodicBoundary)[1][1] == edge[1]
                cable1 = filter(cable->filter(t->edge[1] in t, Weave.periodicBoundary)[1][2] in cable, Weave.cables)[1]
            else
                cable1 = filter(cable->filter(t->edge[1] in t, Weave.periodicBoundary)[1][1] in cable, Weave.cables)[1]
            end
        else
            cable1 = cable1[1]
        end
        if isempty(cable2)
            if filter(t->edge[2] in t, Weave.periodicBoundary)[1][1] == edge[2]
                cable2 = filter(cable->filter(t->edge[2] in t, Weave.periodicBoundary)[1][2] in cable, Weave.cables)[1]
            else
                cable2 = filter(cable->filter(t->edge[2] in t, Weave.periodicBoundary)[1][1] in cable, Weave.cables)[1]
            end
        else
            cable2 = cable2[1]
        end

        v1, v2 = p0[:,edge[2]]-p0[:,edge[1]] + p0[:,cable1[2]]-p0[:,cable1[1]], p0[:,cable2[2]]-p0[:,cable2[1]] + p0[:,edge[2]]-p0[:,edge[1]] 
        n1, n2 = vcat(p0[1:2,edge[1]],0), vcat(p0[1:2,edge[2]],0)
        v1, v2 = v1 - n1 * (n1' * v1), v2 - n2 * (n2' * v2)
        v1 = 0.75*v1; v2= 0.75*v2
        A = [([3 2 1 0; 0 0 1 0; 1 1 1 1; 0 0 0 1] \ [v2[i],v1[i],p0[:,edge[2]][i],p0[:,edge[1]][i]]) for i in 1:3]

        for t in 0:0.05:0.95
            curP = t*(p0[:,edge[2]]-p0[:,edge[1]])+p0[:,edge[1]]
            #push!(newPointList, (outer==1 ? (1+epsilon-(6*epsilon)*t^2+(4*epsilon)*t^3) : (1+epsilon-(6*epsilon)*(1-t)^2+(4*epsilon)*(1-t)^3))*curP./norm(curP))
            push!(newPointList, [A[i]'*[t^3,t^2,t,1] for i in 1:length(A)])
            #push!(newPointList, curP ./ norm(curP))
            push!(newEdgeList, (length(newPointList), length(newPointList)+1))
        end
        push!(newPointList, p0[:,edge[2]])#(1-outer*epsilon)*curP)
    end
    for i in size(p0)[2]
        push!(newPointList, p0[:,i])# ./ norm(p0[:,i]))
    end
    open("weavingframeworkmodeloncylinderkagome.poly", "w") do f
        write(f, "POINTS\n")
        foreach(i->write(f, string("$(i): ", newPointList[i][1], " ", newPointList[i][2], " ", newPointList[i][3],"\n")), 1:length(newPointList))
        write(f,"POLYS\n")
        foreach(i->write(f, string("$(i): ", newEdgeList[i][1], " ", newEdgeList[i][2],"\n")), 1:length(newEdgeList))
        write(f,"END")
    end
end



function test_sphere_tetrahedron()
    Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true, true,false,true,false, false,true,false,true], [(1,5),(2,11),(3,7),(4,9),(6,10),(8,12)], [(1,2),(2,3),(3,4),(4,1), (5,6),(6,7),(7,8),(8,5), (9,10),(10,11),(11,12),(12,9)], []#=[(1,2,3,4), (5,6,7,8), (9,10,11,12)]=#)
    p0 = [1 0 0; 0 -1 0; -1 0 0; 0 1 0; 1 0 0; 0 0 -1; -1 0 0; 0 0 1; 0 1 0; 0 0 -1; 0 -1 0; 0 0 1]'
    initialConfiguration = toArray(p0, Weave) #+ (randn(Float64, length(toArray(p0, Weave))) .- 0.5)*0.05
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-12)
    computeOptimalWeaving(q, Weave)
    normalspace=evaluate.(differentiate(Weave.constraints, Weave.coordinateVariables), Weave.coordinateVariables=>q)
    display("Tangent direction")
    display(toMatrix(collect(nullspace(normalspace)[1:end,1]), Weave))
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:green,:blue,:red])
    toSphere(q, Weave)
end

function test_sphere_rhombicuboctahedron()
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

function test_sphere_a2double()
    Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true, true,false,true,false, false,true,false,true], [(1,5),(2,11),(3,7),(4,9),(6,10),(8,12)], [(1,2),(2,3),(3,4),(4,1), (5,6),(6,7),(7,8),(8,5), (9,10),(10,11),(11,12),(12,9)], []#=[(1,2,3,4), (5,6,7,8), (9,10,11,12)]=#)
    p0 = [1 0 0; 0 -1 0; -1 0 0; 0 1 0; 1 0 0; 0 0 -1; -1 0 0; 0 0 1; 0 1 0; 0 0 -1; 0 -1 0; 0 0 1]'
    p, Weave = doubleRings(p0, Weave)
    initialConfiguration = toArray(p, Weave) #+ (randn(Float64, length(toArray(p0, Weave))) .- 0.5)*0.05
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-12)
    plotWeaving(initialConfiguration, Weave; colorscheme=[:yellow,:cyan,:purple,:yellow,:cyan,:purple])
    q = computeOptimalWeaving(q, Weave; tol=1e-2)
    normalspace=evaluate.(differentiate(Weave.constraints, Weave.coordinateVariables), Weave.coordinateVariables=>q)
    display("Tangent direction")
    display(toMatrix(collect(nullspace(normalspace)[1:end,1]), Weave))
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:yellow,:cyan,:purple])
    toSphere(q, Weave)
end


function test_sphere_octahedron()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false, false,true,false,true,false,true, false,true,false,true,false,true, true,false,true,false,true,false], [(1,7), (4,10), (14,20), (17,23), (3,15), (6,18), (9,16), (12,13), (11,19), (5,24), (2,21), (8,22)], [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1), (7,8),(8,9),(9,10),(10,11),(11,12),(12,7), (13,14),(14,15),(15,16),(16,17),(17,18),(18,13), (19,20),(20,21),(21,22),(22,23),(23,24),(24,19)], []#=[(1,2,3,4,5,6), (7,8,9,10,11,12), (13,14,15,16,17,18), (19,20,21,22,23,24)]=#)
    p0 = [1 0 0; cos(2*pi/6) sin(2*pi/6) 0; cos(4*pi/6) sin(4*pi/6) 0; cos(6*pi/6) sin(6*pi/6) 0; cos(8*pi/6) sin(8*pi/6) 0; cos(10*pi/6) sin(10*pi/6) 0;
        1 0 0; cos(2*pi/6) 0 sin(2*pi/6); cos(4*pi/6) 0 sin(4*pi/6); cos(6*pi/6) 0 sin(6*pi/6); cos(8*pi/6) 0 sin(8*pi/6); cos(10*pi/6) 0 sin(10*pi/6);
        cos(10*pi/6) 0 sin(10*pi/6); 0 sqrt(1/2) -sqrt(1/2); cos(4*pi/6) sin(4*pi/6) 0; cos(4*pi/6) 0 sin(4*pi/6); 0 -sqrt(1/2) sqrt(1/2); cos(10*pi/6) sin(10*pi/6) 0;
        cos(8*pi/6) 0 sin(8*pi/6); 0 sqrt(1/2) -sqrt(1/2); cos(2*pi/6) sin(2*pi/6) 0; cos(2*pi/6) 0 sin(2*pi/6); 0 -sqrt(1/2) sqrt(1/2); cos(8*pi/6) sin(8*pi/6) 0]'
    initialConfiguration = toArray(p0,Weave)

    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave; colorscheme=[:yellow,  :cyan, :green2, :purple2])
    toSphere(q, Weave)
end


function test_spheroid_tetrahedron()
    Q = [-0.9188824523658923 -0.3509183022611896 0.18030913418950417; 0.23167305673407199 -0.11000208491682795 0.9665542592619633; -0.3193471990131362 0.9299225163700875 0.18237730130251661]
    Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true, true,false,true,false, false,true,false,true], [(1,5),(2,11),(3,7),(4,9),(6,10),(8,12)], [(1,2),(2,3),(3,4),(4,1), (5,6),(6,7),(7,8),(8,5), (9,10),(10,11),(11,12),(12,9)], []; rotationsmatrix=Q, manifold="ellipsoid", ellipse=[1,1]#=[(1,2,3,4), (5,6,7,8), (9,10,11,12)]=#)
    p0 = Q*[1 0 0; 0 -1 0; -1 0 0; 0 1 0; 1 0 0; 0 0 -1; -1 0 0; 0 0 1; 0 1 0; 0 0 -1; 0 -1 0; 0 0 1]'
    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    q = computeOptimalWeaving(q, Weave; bending=true)
    #plotWeaving(q, Weave; colorscheme=[:yellow,  :cyan, :green2, :purple2])
    global itercount = 1
    toEllipsoid(q, Weave; ellipse=[1,1], pltstring = "spheroid_tet2_$(itercount)")
    EnergyArray = [(energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; bending=true, no_cables=true, no_struts=true))]
    stepsize = 0.01
    endpoint = 2.0
    for a in vcat(1+stepsize:stepsize:endpoint, endpoint-stepsize:-stepsize:1)
        display(a)
        global itercount = itercount+1
        Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true, true,false,true,false, false,true,false,true], [(1,5),(2,11),(3,7),(4,9),(6,10),(8,12)], [(1,2),(2,3),(3,4),(4,1), (5,6),(6,7),(7,8),(8,5), (9,10),(10,11),(11,12),(12,9)], []; rotationsmatrix=Q, manifold="ellipsoid", ellipse=[a,a]#=[(1,2,3,4), (5,6,7,8), (9,10,11,12)]=#)
        q = newtonCorrect(q, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
        q = computeOptimalWeaving(q, Weave; bending=true)
        #plotWeaving(q, Weave; colorscheme=[:yellow,  :cyan, :green2, :purple2])
        toEllipsoid(q, Weave; ellipse=[a,a], pltstring = "spheroid_tet2_$(itercount)")
        push!(EnergyArray, (energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; bending=true, no_cables=true, no_struts=true)))
        write_txt_file(EnergyArray; filename="Tetrahedron3")
    end
    plot_energy(EnergyArray, Weave; pltlabel="Tetrahedron2", endpoint=endpoint, stepsize=stepsize, bending=true)
    write_txt_file(EnergyArray; filename="Tetrahedron3")
end


function test_spheroid_octahedron()
    Q = [-0.9188824523658923 -0.3509183022611896 0.18030913418950417; 0.23167305673407199 -0.11000208491682795 0.9665542592619633; -0.3193471990131362 0.9299225163700875 0.18237730130251661]
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false, false,true,false,true,false,true, false,true,false,true,false,true, true,false,true,false,true,false], [(1,7), (4,10), (14,20), (17,23), (3,15), (6,18), (9,16), (12,13), (11,19), (5,24), (2,21), (8,22)], [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1), (7,8),(8,9),(9,10),(10,11),(11,12),(12,7), (13,14),(14,15),(15,16),(16,17),(17,18),(18,13), (19,20),(20,21),(21,22),(22,23),(23,24),(24,19)], []; rotationsmatrix=Q, manifold="ellipsoid", ellipse=[1,1]#=[(1,2,3,4,5,6), (7,8,9,10,11,12), (13,14,15,16,17,18), (19,20,21,22,23,24)]=#)
    p0 = Q*[1 0 0; cos(2*pi/6) sin(2*pi/6) 0; cos(4*pi/6) sin(4*pi/6) 0; cos(6*pi/6) sin(6*pi/6) 0; cos(8*pi/6) sin(8*pi/6) 0; cos(10*pi/6) sin(10*pi/6) 0;
        1 0 0; cos(2*pi/6) 0 sin(2*pi/6); cos(4*pi/6) 0 sin(4*pi/6); cos(6*pi/6) 0 sin(6*pi/6); cos(8*pi/6) 0 sin(8*pi/6); cos(10*pi/6) 0 sin(10*pi/6);
        cos(10*pi/6) 0 sin(10*pi/6); 0 sqrt(1/2) -sqrt(1/2); cos(4*pi/6) sin(4*pi/6) 0; cos(4*pi/6) 0 sin(4*pi/6); 0 -sqrt(1/2) sqrt(1/2); cos(10*pi/6) sin(10*pi/6) 0;
        cos(8*pi/6) 0 sin(8*pi/6); 0 sqrt(1/2) -sqrt(1/2); cos(2*pi/6) sin(2*pi/6) 0; cos(2*pi/6) 0 sin(2*pi/6); 0 -sqrt(1/2) sqrt(1/2); cos(8*pi/6) sin(8*pi/6) 0]'
    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    q = computeOptimalWeaving(q, Weave; bending = true)
    #plotWeaving(q, Weave; colorscheme=[:yellow,  :cyan, :green2, :purple2])
    global itercount = 1
    toEllipsoid(q, Weave; ellipse=[1,1], pltstring = "spheroid_oct_$(itercount)")
    EnergyArray = [(energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true))]
    endpoint = 2
    stepsize = 0.01
    for a in vcat(1+stepsize:stepsize:endpoint, endpoint-stepsize:-stepsize:1)
        display(a)
        global itercount = itercount+1
        Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false, false,true,false,true,false,true, false,true,false,true,false,true, true,false,true,false,true,false], [(1,7), (4,10), (14,20), (17,23), (3,15), (6,18), (9,16), (12,13), (11,19), (5,24), (2,21), (8,22)], [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1), (7,8),(8,9),(9,10),(10,11),(11,12),(12,7), (13,14),(14,15),(15,16),(16,17),(17,18),(18,13), (19,20),(20,21),(21,22),(22,23),(23,24),(24,19)], []; rotationsmatrix=Q, manifold="ellipsoid", ellipse=[a,a]#=[(1,2,3,4,5,6), (7,8,9,10,11,12), (13,14,15,16,17,18), (19,20,21,22,23,24)]=#)
        q = newtonCorrect(q, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
        q = computeOptimalWeaving(q, Weave; bending = true)
        #plotWeaving(q, Weave; colorscheme=[:yellow,  :cyan, :green2, :purple2])
        toEllipsoid(q, Weave; ellipse=[a,a], pltstring = "spheroid_oct_$(itercount)")
        push!(EnergyArray, (energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true)))
    end
    plot_energy(EnergyArray, Weave; pltlabel="Octahedron", bending=true)
    write_txt_file(EnergyArray; filename="Octahedron3")
end

function test_spheroid_rhombioctahedron()
    Q = [-0.9188824523658923 -0.3509183022611896 0.18030913418950417; 0.23167305673407199 -0.11000208491682795 0.9665542592619633; -0.3193471990131362 0.9299225163700875 0.18237730130251661]
    Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true,false,true,false,true, true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true, true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true], 
            [(4, 43), (9, 28), (15, 47), (29, 38), (8, 40), (2, 19), (13, 23), (12, 44), (21, 37), (26, 33), (30, 46), (22, 45), (11, 36), (7, 48), (3, 35), (14, 31), (17, 42), (6, 32), (16, 39), (5, 24), (18, 34), (10, 20), (25, 41), (1, 27)], 
            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,1), (9,10),(10,11),(11,12),(12,13),(13,14),(14,15),(15,16),(16,9), (17,18),(18,19),(19,20),(20,21),(21,22),(22,23),(23,24),(24,17), (25,26),(26,27),(27,28),(28,29),(29,30),(30,31),(31,32),(32,25), (33,34),(34,35),(35,36),(36,37),(37,38),(38,39),(39,40),(40,33), (41,42),(42,43),(43,44),(44,45),(45,46),(46,47),(47,48),(48,41)], []; 
            rotationsmatrix=Q, offset = 0.1, manifold="ellipsoid", ellipse=[1,1])
    p0 = Q*[-1/3 sqrt(7/9) 1/3; 1/3 sqrt(7/9) 1/3; sqrt(7/9) 1/3 1/3; sqrt(7/9) -1/3 1/3; 1/3 -sqrt(7/9) 1/3; -1/3 -sqrt(7/9) 1/3; -sqrt(7/9) -1/3 1/3; -sqrt(7/9) 1/3 1/3;
    -1/3 sqrt(7/9) -1/3; 1/3 sqrt(7/9) -1/3; sqrt(7/9) 1/3 -1/3; sqrt(7/9) -1/3 -1/3; 1/3 -sqrt(7/9) -1/3; -1/3 -sqrt(7/9) -1/3; -sqrt(7/9) -1/3 -1/3; -sqrt(7/9) 1/3 -1/3;
    1/3 -1/3 sqrt(7/9); 1/3 1/3 sqrt(7/9); 1/3 sqrt(7/9) 1/3; 1/3 sqrt(7/9) -1/3; 1/3 1/3 -sqrt(7/9); 1/3 -1/3 -sqrt(7/9); 1/3 -sqrt(7/9) -1/3; 1/3 -sqrt(7/9) 1/3;
    -1/3 -1/3 sqrt(7/9); -1/3 1/3 sqrt(7/9); -1/3 sqrt(7/9) 1/3; -1/3 sqrt(7/9) -1/3; -1/3 1/3 -sqrt(7/9); -1/3 -1/3 -sqrt(7/9); -1/3 -sqrt(7/9) -1/3; -1/3 -sqrt(7/9) 1/3;
    -1/3 1/3 sqrt(7/9); 1/3 1/3 sqrt(7/9); sqrt(7/9) 1/3 1/3; sqrt(7/9) 1/3 -1/3; 1/3 1/3 -sqrt(7/9); -1/3 1/3 -sqrt(7/9); -sqrt(7/9) 1/3 -1/3; -sqrt(7/9) 1/3 1/3;
    -1/3 -1/3 sqrt(7/9); 1/3 -1/3 sqrt(7/9); sqrt(7/9) -1/3 1/3; sqrt(7/9) -1/3 -1/3; 1/3 -1/3 -sqrt(7/9); -1/3 -1/3 -sqrt(7/9); -sqrt(7/9) -1/3 -1/3; -sqrt(7/9) -1/3 1/3;
    ]'
    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    q = computeOptimalWeaving(q, Weave; bending = false)
    #plotWeaving(q, Weave; colorscheme=[:yellow,  :cyan, :green2, :purple2])
    global itercount = 1
    toEllipsoid(q, Weave; ellipse=[1,1], pltstring = "spheroid_rhomboct_$(itercount)")
    EnergyArray = [(energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true))]
    endpoint = 2
    stepsize = 0.01
    for a in vcat(1+stepsize:stepsize:endpoint, endpoint-stepsize:-stepsize:1)
        display(a)
        global itercount = itercount+1
        Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true,false,true,false,true, true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true, true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true], 
            [(4, 43), (9, 28), (15, 47), (29, 38), (8, 40), (2, 19), (13, 23), (12, 44), (21, 37), (26, 33), (30, 46), (22, 45), (11, 36), (7, 48), (3, 35), (14, 31), (17, 42), (6, 32), (16, 39), (5, 24), (18, 34), (10, 20), (25, 41), (1, 27)], 
            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,1), (9,10),(10,11),(11,12),(12,13),(13,14),(14,15),(15,16),(16,9), (17,18),(18,19),(19,20),(20,21),(21,22),(22,23),(23,24),(24,17), (25,26),(26,27),(27,28),(28,29),(29,30),(30,31),(31,32),(32,25), (33,34),(34,35),(35,36),(36,37),(37,38),(38,39),(39,40),(40,33), (41,42),(42,43),(43,44),(44,45),(45,46),(46,47),(47,48),(48,41)], []; 
            rotationsmatrix=Q, offset = 0.1, manifold="ellipsoid", ellipse=[a,a])
        q = newtonCorrect(q, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
        q = computeOptimalWeaving(q, Weave; bending = false)
        #plotWeaving(q, Weave; colorscheme=[:yellow,  :cyan, :green2, :purple2])
        toEllipsoid(q, Weave; ellipse=[a,a], pltstring = "spheroid_rhomboct_$(itercount)")
        push!(EnergyArray, (energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true)))
        write_txt_file(EnergyArray; filename="RhombiOctahedron3")
    end
    #plot_energy(EnergyArray, Weave; pltlabel="RhombiOctahedron", bending=false)
end

function test_spheroid_dodecahedron()
    Q = [-0.9188824523658923 -0.3509183022611896 0.18030913418950417; 0.23167305673407199 -0.11000208491682795 0.9665542592619633; -0.3193471990131362 0.9299225163700875 0.18237730130251661]
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true,false,true, false,true,false,true,false,true,false,true,false,true, true,false,true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true,false,true], [(5, 59), (14, 23), (19, 28), (25, 47), (34, 45), (41, 51), (2, 43), (35, 55), (24, 57), (3, 21), (13, 33), (12, 44), (6, 16), (39, 50), (27, 37), (29, 52), (1, 11), (18, 38), (7, 48), (20, 53), (9, 36), (17, 49), (40, 60), (30, 42), (10, 54), (22, 32), (15, 58), (46, 56), (8, 26), (4, 31)], 
                                            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(10,1), (11,12),(12,13),(13,14),(14,15),(15,16),(16,17),(17,18),(18,19),(19,20),(20,11), (21,22),(22,23),(23,24),(24,25),(25,26),(26,27),(27,28),(28,29),(29,30),(30,21), (31,32),(32,33),(33,34),(34,35),(35,36),(36,37),(37,38),(38,39),(39,40),(40,31), (41,42),(42,43),(43,44),(44,45),(45,46),(46,47),(47,48),(48,49),(49,50),(50,41), (51,52),(52,53),(53,54),(54,55),(55,56),(56,57),(57,58),(58,59),(59,60),(60,51)], 
                                            []#=[(1,2,3,4,5,6,7,8,9,10), (11,12,13,14,15,16,17,18,19,20), (21,22,23,24,25,26,27,28,29,30), (31,32,33,34,35,36,37,38,39,40), (41,42,43,44,45,46,47,48,49,50), (51,52,53,54,55,56,57,58,59,60)]=#; 
                                            rotationsmatrix=Q, manifold="ellipsoid", ellipse=[1,1], offset=0.1)
    p0 = Q*[1 0 0; cos(2*pi/10) sin(2*pi/10) 0; cos(4*pi/10) sin(4*pi/10) 0; cos(6*pi/10) sin(6*pi/10) 0; cos(8*pi/10) sin(8*pi/10) 0; cos(10*pi/10) sin(10*pi/10) 0; cos(12*pi/10) sin(12*pi/10) 0; cos(14*pi/10) sin(14*pi/10) 0; cos(16*pi/10) sin(16*pi/10) 0; cos(18*pi/10) sin(18*pi/10) 0;
          1 0 0; cos(2*pi/10) 0 sin(2*pi/10); cos(4*pi/10) 0 sin(4*pi/10); cos(6*pi/10) 0 sin(6*pi/10); cos(8*pi/10) 0 sin(8*pi/10); cos(10*pi/10) 0 sin(10*pi/10); cos(12*pi/10) 0 sin(12*pi/10); cos(14*pi/10) 0 sin(14*pi/10); cos(16*pi/10) 0 sin(16*pi/10); cos(18*pi/10) 0 sin(18*pi/10);
          cos(4*pi/10) sin(4*pi/10) 0; 0 sqrt(1/2) sqrt(1/2); cos(6*pi/10) 0 sin(6*pi/10); -1/2 -sqrt(1/8) sqrt(5/8); -0.4 -sqrt(0.5) sqrt(0.34); cos(14*pi/10) sin(14*pi/10) 0; 0 -sqrt(1/2) -sqrt(1/2); cos(16*pi/10) 0 sin(16*pi/10); 1/2 sqrt(1/8) -sqrt(5/8); 1/2 sqrt(3/8) -sqrt(3/8);
          cos(6*pi/10) sin(6*pi/10) 0; 0 sqrt(1/2) sqrt(1/2); cos(4*pi/10) 0 sin(4*pi/10); 1/2 -sqrt(1/8) sqrt(5/8); 0.4 -sqrt(0.5) sqrt(0.34); cos(16*pi/10) sin(16*pi/10) 0; 0 -sqrt(1/2) -sqrt(1/2); cos(14*pi/10) 0 sin(14*pi/10); -1/2 sqrt(1/8) -sqrt(5/8); -1/2 sqrt(3/8) -sqrt(3/8);
          0 sqrt(3/8) -sqrt(5/8); 1/2 sqrt(3/8) -sqrt(3/8); cos(2*pi/10) sin(2*pi/10) 0; cos(2*pi/10) 0 sin(2*pi/10); 1/2 -sqrt(1/8) sqrt(5/8); 0 -sqrt(3/8) sqrt(5/8); -0.4 -sqrt(0.5) sqrt(0.34); cos(12*pi/10) sin(12*pi/10) 0; cos(12*pi/10) 0 sin(12*pi/10); -1/2 sqrt(1/8) -sqrt(5/8);
          0 sqrt(3/8) -sqrt(5/8); 1/2 sqrt(1/8) -sqrt(5/8); cos(18*pi/10) 0 sin(18*pi/10); cos(18*pi/10) sin(18*pi/10) 0; 0.4 -sqrt(0.5) sqrt(0.34); 0 -sqrt(3/8) sqrt(5/8); -1/2 -sqrt(1/8) sqrt(5/8); cos(8*pi/10) 0 sin(8*pi/10); cos(8*pi/10) sin(8*pi/10) 0; -1/2 sqrt(3/8) -sqrt(3/8)
    ]'

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    q = computeOptimalWeaving(q, Weave; bending=false)
    global itercount = 1
    #plotWeaving(q, Weave; colorscheme=[:yellow,  :cyan, :green2, :purple2, :magenta2, :blue2])
    toEllipsoid(q, Weave; ellipse=[1,1], pltstring = "spheroid_dod_$(itercount)")
    EnergyArray = [(energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true))]
    endpoint = 2
    stepsize = 0.01
    for a in vcat(1+stepsize:stepsize:endpoint, endpoint-stepsize:-stepsize:1)
        display(a)
        global itercount = itercount+1
        Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true,false,true, false,true,false,true,false,true,false,true,false,true, true,false,true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true,false,true], [(5, 59), (14, 23), (19, 28), (25, 47), (34, 45), (41, 51), (2, 43), (35, 55), (24, 57), (3, 21), (13, 33), (12, 44), (6, 16), (39, 50), (27, 37), (29, 52), (1, 11), (18, 38), (7, 48), (20, 53), (9, 36), (17, 49), (40, 60), (30, 42), (10, 54), (22, 32), (15, 58), (46, 56), (8, 26), (4, 31)], 
                                                        [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(10,1), (11,12),(12,13),(13,14),(14,15),(15,16),(16,17),(17,18),(18,19),(19,20),(20,11), (21,22),(22,23),(23,24),(24,25),(25,26),(26,27),(27,28),(28,29),(29,30),(30,21), (31,32),(32,33),(33,34),(34,35),(35,36),(36,37),(37,38),(38,39),(39,40),(40,31), (41,42),(42,43),(43,44),(44,45),(45,46),(46,47),(47,48),(48,49),(49,50),(50,41), (51,52),(52,53),(53,54),(54,55),(55,56),(56,57),(57,58),(58,59),(59,60),(60,51)], 
                                                        []#=[(1,2,3,4,5,6,7,8,9,10), (11,12,13,14,15,16,17,18,19,20), (21,22,23,24,25,26,27,28,29,30), (31,32,33,34,35,36,37,38,39,40), (41,42,43,44,45,46,47,48,49,50), (51,52,53,54,55,56,57,58,59,60)]=#; 
                                                        rotationsmatrix=Q, manifold="ellipsoid", ellipse=[a,a], offset=0.1)
        q = newtonCorrect(q, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
        q = computeOptimalWeaving(q, Weave; bending=false)
        #plotWeaving(q, Weave; colorscheme=[:yellow,  :cyan, :green2, :purple2, :magenta2, :blue2])
        toEllipsoid(q, Weave; ellipse=[a,a], pltstring = "spheroid_dod_$(itercount)")
        push!(EnergyArray, (energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true)))
        write_txt_file(EnergyArray; filename="Dodecahedron3")
    end
    #plot_energy(EnergyArray, Weave; pltlabel="Dodecahedron", bending=false)
end

function test_spheroid_trunctetra()
    Q = [-0.9188824523658923 -0.3509183022611896 0.18030913418950417; 0.23167305673407199 -0.11000208491682795 0.9665542592619633; -0.3193471990131362 0.9299225163700875 0.18237730130251661]
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false, false,true,false,true,false,true, false,true,false,true,false,true, false,true,false,true,false,true], [(1,13),(2,8),(3,9),(4,24),(5,19),(6,18),(7,14),(10,23),(11,22),(12,15),(16,21),(17,20)], 
                                            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1), (7,8),(8,9),(9,10),(10,11),(11,12),(12,7), (13,14),(14,15),(15,16),(16,17),(17,18),(18,13), (19,20),(20,21),(21,22),(22,23),(23,24),(24,19)], 
                                            []#=[(1,2,3,4,5,6,7,8,9,10), (11,12,13,14,15,16,17,18,19,20), (21,22,23,24,25,26,27,28,29,30), (31,32,33,34,35,36,37,38,39,40), (41,42,43,44,45,46,47,48,49,50), (51,52,53,54,55,56,57,58,59,60)]=#; 
                                            rotationsmatrix=Q, ellipse=[1,1], manifold = "ellipsoid")
    p0 = Q*[1 3 1; 1 1 3; -1 -1 3; -3 -1 1; -3 1 -1; -1 3 -1; 
          3 1 1; 1 1 3; -1 -1 3; -1 -3 1; 1 -3 -1; 3 -1 -1; 
          1 3 1; 3 1 1; 3 -1 -1; 1 -1 -3; -1 1 -3; -1 3 -1;
          -3 1 -1; -1 1 -3; 1 -1 -3; 1 -3 -1; -1 -3 1; -3 -1 1;
    ]' ./ sqrt(11)

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    q = computeOptimalWeaving(q, Weave; bending=true)
    global itercount = 1
    toEllipsoid(q, Weave; ellipse=[1,1], pltstring = "spheroid_trunctet_$(itercount)")
    EnergyArray = [(energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true))]
    endpoint = 2
    stepsize = 0.01
    for a in vcat(1+stepsize:stepsize:endpoint, endpoint-stepsize:-stepsize:1)
        display(a)
        global itercount = itercount+1
        Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false, false,true,false,true,false,true, false,true,false,true,false,true, false,true,false,true,false,true], [(1,13),(2,8),(3,9),(4,24),(5,19),(6,18),(7,14),(10,23),(11,22),(12,15),(16,21),(17,20)], 
            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1), (7,8),(8,9),(9,10),(10,11),(11,12),(12,7), (13,14),(14,15),(15,16),(16,17),(17,18),(18,13), (19,20),(20,21),(21,22),(22,23),(23,24),(24,19)], 
            []#=[(1,2,3,4,5,6,7,8,9,10), (11,12,13,14,15,16,17,18,19,20), (21,22,23,24,25,26,27,28,29,30), (31,32,33,34,35,36,37,38,39,40), (41,42,43,44,45,46,47,48,49,50), (51,52,53,54,55,56,57,58,59,60)]=#; 
            rotationsmatrix=Q, ellipse=[a,a], manifold = "ellipsoid")
        q = newtonCorrect(q, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
        q = computeOptimalWeaving(q, Weave; bending=true)
        toEllipsoid(q, Weave; ellipse=[a,a], pltstring = "spheroid_trunctet_$(itercount)")
        push!(EnergyArray, (energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true)))
        write_txt_file(EnergyArray; filename="TruncTetra3")
    end
    #plot_energy(EnergyArray, Weave; pltlabel="TruncTetra", bending=false)
end

function test_spheroid_truncatedcube()
    offset = 0.1
    Q = [-0.9188824523658923 -0.3509183022611896 0.18030913418950417; 0.23167305673407199 -0.11000208491682795 0.9665542592619633; -0.3193471990131362 0.9299225163700875 0.18237730130251661]
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true], [(1, 40), (2, 39), (3, 16), (4, 15), (5, 43), (6, 42), (7, 28), (8, 27), (9, 38), (10, 37), (11, 24), (12, 23), (13, 45), (14, 44), (17, 36), (18, 35), (19, 32), (20, 31), (21, 47), (22, 46), (25, 34), (26, 33), (29, 41), (30, 48)], 
                                            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,1), (9,10),(10,11),(11,12),(12,13),(13,14),(14,15),(15,16),(16,9), (17,18),(18,19),(19,20),(20,21),(21,22),(22,23),(23,24),(24,17), (25,26),(26,27),(27,28),(28,29),(29,30),(30,31),(31,32),(32,25), (33,34),(34,35),(35,36),(36,37),(37,38),(38,39),(39,40),(40,33), (41,42),(42,43),(43,44),(44,45),(45,46),(46,47),(47,48),(48,41)], 
                                            []#=[(1,2,3,4,5,6,7,8,9,10), (11,12,13,14,15,16,17,18,19,20), (21,22,23,24,25,26,27,28,29,30), (31,32,33,34,35,36,37,38,39,40), (41,42,43,44,45,46,47,48,49,50), (51,52,53,54,55,56,57,58,59,60)]=#; 
                                            rotationsmatrix=Q, manifold="ellipsoid", samples = [Q*[sqrt(2)-1,1,1]])
    ξ = sqrt(2)-1
    p0 = Q*[ξ 1 1; -ξ 1 1; -1 ξ 1; -1 -ξ 1; -ξ -1 1; ξ -1 1; 1 -ξ 1; 1 ξ 1;
        -1 1 ξ; -1 1 -ξ; -1 ξ -1; -1 -ξ -1; -1 -1 -ξ; -1 -1 ξ; -1 -ξ 1; -1 ξ 1;
        -ξ 1 -1; ξ 1 -1; 1 ξ -1; 1 -ξ -1; ξ -1 -1; -ξ -1 -1; -1 -ξ -1; -1 ξ -1;
        1 1 -ξ; 1 1 ξ; 1 ξ 1; 1 -ξ 1; 1 -1 ξ; 1 -1 -ξ; 1 -ξ -1; 1 ξ -1;
        1 1 ξ; 1 1 -ξ; ξ 1 -1; -ξ 1 -1; -1 1 -ξ; -1 1 ξ; -ξ 1 1; ξ 1 1;
        1 -1 ξ; ξ -1 1; -ξ -1 1; -1 -1 ξ; -1 -1 -ξ; -ξ -1 -1; ξ -1 -1; 1 -1 -ξ;
    ]' ./ sqrt(ξ^2+2)

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-10)
    q = computeOptimalWeaving(q, Weave; bending=false)
    global itercount = 1
    toEllipsoid(q, Weave; ellipse=[1,1], pltstring = "spheroid_trunccub_$(itercount)")
    EnergyArray = [(energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true))]
    endpoint = 2
    stepsize = 0.01
    for a in vcat(1+stepsize:stepsize:endpoint, endpoint-stepsize:-stepsize:1)
        display(a)
        global itercount = itercount+1
        Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true], [(1, 40), (2, 39), (3, 16), (4, 15), (5, 43), (6, 42), (7, 28), (8, 27), (9, 38), (10, 37), (11, 24), (12, 23), (13, 45), (14, 44), (17, 36), (18, 35), (19, 32), (20, 31), (21, 47), (22, 46), (25, 34), (26, 33), (29, 41), (30, 48)], 
            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,1), (9,10),(10,11),(11,12),(12,13),(13,14),(14,15),(15,16),(16,9), (17,18),(18,19),(19,20),(20,21),(21,22),(22,23),(23,24),(24,17), (25,26),(26,27),(27,28),(28,29),(29,30),(30,31),(31,32),(32,25), (33,34),(34,35),(35,36),(36,37),(37,38),(38,39),(39,40),(40,33), (41,42),(42,43),(43,44),(44,45),(45,46),(46,47),(47,48),(48,41)], 
            []#=[(1,2,3,4,5,6,7,8,9,10), (11,12,13,14,15,16,17,18,19,20), (21,22,23,24,25,26,27,28,29,30), (31,32,33,34,35,36,37,38,39,40), (41,42,43,44,45,46,47,48,49,50), (51,52,53,54,55,56,57,58,59,60)]=#; 
            rotationsmatrix=Q, ellipse=[a,a], manifold="ellipsoid", samples = [Q*[sqrt(2)-1,1,1]])
        q = newtonCorrect(q, Weave.coordinateVariables, Weave.constraints; tol = 1e-10)
        q = computeOptimalWeaving(q, Weave; bending=false)
        toEllipsoid(q, Weave; ellipse=[a,a], pltstring = "spheroid_trunccub_$(itercount)")
        push!(EnergyArray, (energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true)))
        write_txt_file(EnergyArray; filename="TruncCube3")
    end
    #plot_energy(EnergyArray, Weave; pltlabel="TruncCube", bending=false)
end

function test_spheroid_truncatedoctahedron()
    offset = 0.1
    Q = [-0.9188824523658923 -0.3509183022611896 0.18030913418950417; 0.23167305673407199 -0.11000208491682795 0.9665542592619633; -0.3193471990131362 0.9299225163700875 0.18237730130251661]
    pointlist, edgelist = read_poly_file("trunc_oct.poly")
    newedgelist = [edgelist[1]]
    deleteat!(edgelist, 1)
    while ! isempty(edgelist)
        index = findfirst(ring -> any(edge -> (edge in ring)||((edge[2],edge[1]) in ring), vcat(newedgelist...)), edgelist)
        push!(newedgelist, popat!(edgelist, index))
    end
    cables = []
    for ring in newedgelist
        ring_sugg = [(ring[i][1]+mod(i,2)*length(pointlist), ring[i][2]+mod(i+1,2)*length(pointlist)) for i in 1:length(ring)]
        if !any(ringlist -> any(cable->ringlist==cable||ringlist==(cable[2],cable[1]), cables), ring_sugg)
            #println(filter(t->any(entry->entry in t, vcat([[r[1], r[2]] for r in ring_sugg]...)), cables), ring_sugg)
            append!(cables, ring_sugg)
        else
            ring_sugg = [(ring[i][1]+mod(i+1,2)*length(pointlist), ring[i][2]+mod(i,2)*length(pointlist)) for i in 1:length(ring)]
            if any(ringlist -> any(cable->ringlist==cable||ringlist==(cable[2],cable[1]), cables), ring_sugg)
                continue
            end    
            append!(cables, ring_sugg)
        end
        println("")
    end
    p0 = Q*hcat(vcat(pointlist .* (1+offset), pointlist .* (1-offset))...)
    barlist = [(bar, bar+length(pointlist)) for bar in 1:length(pointlist)]
    offsetlist = vcat([true for _ in 1:length(pointlist)], [false for _ in 1:length(pointlist)])
    Weave = WeavingsOnManifolds.WeavingOnManifold(offsetlist, barlist, cables, []; 
        rotationsmatrix=Q, samples = [Q*pointlist[1]], manifold="ellipsoid", ellipse=[1,1])

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-10)
    q = computeOptimalWeaving(q, Weave; bending=false)
    global itercount = 1
    toEllipsoid(q, Weave; ellipse=[1,1], pltstring = "spheroid_truncoct_$(itercount)")
    EnergyArray = [(energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true))]
    endpoint = 2
    stepsize = 0.01
    for a in vcat(1+stepsize:stepsize:endpoint, endpoint-stepsize:-stepsize:1)
        display(a)
        global itercount = itercount+1
        Weave = WeavingsOnManifolds.WeavingOnManifold(offsetlist, barlist, cables, []; 
            rotationsmatrix=Q, samples = [Q*pointlist[1]], manifold="ellipsoid", ellipse=[a,a])
        q = newtonCorrect(q, Weave.coordinateVariables, Weave.constraints; tol = 1e-10)
        q = computeOptimalWeaving(q, Weave; bending=false)
        toEllipsoid(q, Weave; ellipse=[a,a], pltstring = "spheroid_truncoct_$(itercount)")
        push!(EnergyArray, (energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true)))
        write_txt_file(EnergyArray; filename="TruncOct3")
    end
    #plot_energy(EnergyArray, Weave; pltlabel="TruncOct", bending=false)
end

function test_spheroid_truncatedicosihedron()
    offset = 0.1
    Q = [-0.9188824523658923 -0.3509183022611896 0.18030913418950417; 0.23167305673407199 -0.11000208491682795 0.9665542592619633; -0.3193471990131362 0.9299225163700875 0.18237730130251661]
    pointlist, edgelist = read_poly_file("trunc_icos.poly")
    display(pointlist)
    newedgelist = [edgelist[1]]
    deleteat!(edgelist, 1)
    while ! isempty(edgelist)
        index = findfirst(ring -> any(edge -> (edge in ring)||((edge[2],edge[1]) in ring), vcat(newedgelist...)), edgelist)
        push!(newedgelist, popat!(edgelist, index))
    end
    cables = []
    for ring in newedgelist
        ring_sugg = [(ring[i][1]+mod(i,2)*length(pointlist), ring[i][2]+mod(i+1,2)*length(pointlist)) for i in 1:length(ring)]
        if !any(ringlist -> any(cable->ringlist==cable||ringlist==(cable[2],cable[1]), cables), ring_sugg)
            #println(filter(t->any(entry->entry in t, vcat([[r[1], r[2]] for r in ring_sugg]...)), cables), ring_sugg)
            append!(cables, ring_sugg)
        else
            ring_sugg = [(ring[i][1]+mod(i+1,2)*length(pointlist), ring[i][2]+mod(i,2)*length(pointlist)) for i in 1:length(ring)]
            if any(ringlist -> any(cable->ringlist==cable||ringlist==(cable[2],cable[1]), cables), ring_sugg)
                continue
            end    
            append!(cables, ring_sugg)
        end
        println("")
    end
    p0 = Q*hcat(vcat(pointlist .* (1+offset), pointlist .* (1-offset))...)
    barlist = [(bar, bar+length(pointlist)) for bar in 1:length(pointlist)]
    offsetlist = vcat([true for _ in 1:length(pointlist)], [false for _ in 1:length(pointlist)])
    Weave = WeavingsOnManifolds.WeavingOnManifold(offsetlist, barlist, cables, []; 
        rotationsmatrix=Q, manifold = "ellipsoid", ellipse=[1,1], samples = [Q * pointlist[1] ./ norm(pointlist[1])])

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-10)
    q = computeOptimalWeaving(q, Weave; bending=false)
    global itercount = 1
    toEllipsoid(q, Weave; ellipse=[1,1], pltstring = "spheroid_truncico_$(itercount)")
    EnergyArray = [(energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true))]
    write_txt_file(EnergyArray; filename="TruncIcos3")
    endpoint = 2
    stepsize = 0.01
    for a in vcat(1+stepsize:stepsize:endpoint, endpoint-stepsize:-stepsize:1)
        display(a)
        global itercount = itercount+1
        Weave = WeavingsOnManifolds.WeavingOnManifold(offsetlist, barlist, cables, []; 
            rotationsmatrix=Q, manifold = "ellipsoid", ellipse=[a,a], samples = [Q * pointlist[1] ./ norm(pointlist[1])])
        q = newtonCorrect(q, Weave.coordinateVariables, Weave.constraints; tol = 1e-10)
        q = computeOptimalWeaving(q, Weave; bending=false)
        toEllipsoid(q, Weave; ellipse=[a,a], pltstring = "spheroid_truncico_$(itercount)")
        push!(EnergyArray, (energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true)))
        write_txt_file(EnergyArray; filename="TruncIcos3")
    end
    #plot_energy(EnergyArray, Weave; pltlabel="TruncIcos", bending=false)
    #write_txt_file(EnergyArray; filename="TruncIcos")
end

function test_spheroid_truncateddodecahedron()
    Q = [-0.9188824523658923 -0.3509183022611896 0.18030913418950417; 0.23167305673407199 -0.11000208491682795 0.9665542592619633; -0.3193471990131362 0.9299225163700875 0.18237730130251661]
    offset = 0.1
    pointlist, edgelist = read_poly_file("trunc_dodec.poly")
    newedgelist = [edgelist[1]]
    deleteat!(edgelist, 1)
    while ! isempty(edgelist)
        index = findfirst(ring -> any(edge -> (edge in ring)||((edge[2],edge[1]) in ring), vcat(newedgelist...)), edgelist)
        push!(newedgelist, popat!(edgelist, index))
    end
    cables = []
    for ring in newedgelist
        ring_sugg = [(ring[i][1]+mod(i,2)*length(pointlist), ring[i][2]+mod(i+1,2)*length(pointlist)) for i in 1:length(ring)]
        if !any(ringlist -> any(cable->ringlist==cable||ringlist==(cable[2],cable[1]), cables), ring_sugg)
            #println(filter(t->any(entry->entry in t, vcat([[r[1], r[2]] for r in ring_sugg]...)), cables), ring_sugg)
            append!(cables, ring_sugg)
        else
            ring_sugg = [(ring[i][1]+mod(i+1,2)*length(pointlist), ring[i][2]+mod(i,2)*length(pointlist)) for i in 1:length(ring)]
            if any(ringlist -> any(cable->ringlist==cable||ringlist==(cable[2],cable[1]), cables), ring_sugg)
                continue
            end
            append!(cables, ring_sugg)
        end
    end
    p0 = Q*hcat(vcat(pointlist .* (1+offset), pointlist .* (1-offset))...)
    barlist = [(bar, bar+length(pointlist)) for bar in 1:length(pointlist)]
    offsetlist = vcat([true for _ in 1:length(pointlist)], [false for _ in 1:length(pointlist)])

    Weave = WeavingsOnManifolds.WeavingOnManifold(offsetlist, barlist, cables, []; 
            ellipse = [1,1], rotationsmatrix=Q, samples = [Q*pointlist[1]], manifold="ellipsoid")

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    q = computeOptimalWeaving(q, Weave; bending=false)
    global itercount = 1
    toEllipsoid(q, Weave; ellipse=[1,1], pltstring = "spheroid_truncdod_$(itercount)")
    EnergyArray = [(energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true))]
    endpoint = 2
    stepsize = 0.01
    write_txt_file(EnergyArray; filename="TruncDod3")
    for a in vcat(1+stepsize:stepsize:endpoint, endpoint-stepsize:-stepsize:1)
        display(a)
        global itercount = itercount+1
        Weave = WeavingsOnManifolds.WeavingOnManifold(offsetlist, barlist, cables, []; 
                ellipse = [a,a], rotationsmatrix=Q, samples = [Q*pointlist[1]], manifold="ellipsoid")
        q = newtonCorrect(q, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
        q = computeOptimalWeaving(q, Weave; bending=false)
        toEllipsoid(q, Weave; ellipse=[a,a], pltstring = "spheroid_truncdod_$(itercount)")
        push!(EnergyArray, (energyFunction(q, Weave; bending=true), energyFunction(q, Weave; no_struts=true), energyFunction(q, Weave; no_cables=true), energyFunction(q, Weave; no_struts=true, no_cables=true, bending=true)))
        write_txt_file(EnergyArray; filename="TruncDod3")
    end
    #plot_energy(EnergyArray, Weave; pltlabel="TruncDod", bending=false)
    #write_txt_file(EnergyArray; filename="TruncDod")
end










function test_sphere_b2double()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false, false,true,false,true,false,true, false,true,false,true,false,true, true,false,true,false,true,false], [(1,7), (4,10), (14,20), (17,23), (3,15), (6,18), (9,16), (12,13), (11,19), (5,24), (2,21), (8,22)], [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1), (7,8),(8,9),(9,10),(10,11),(11,12),(12,7), (13,14),(14,15),(15,16),(16,17),(17,18),(18,13), (19,20),(20,21),(21,22),(22,23),(23,24),(24,19)], []#=[(1,2,3,4,5,6), (7,8,9,10,11,12), (13,14,15,16,17,18), (19,20,21,22,23,24)]=#)
    p0 = [1 0 0; cos(2*pi/6) sin(2*pi/6) 0; cos(4*pi/6) sin(4*pi/6) 0; cos(6*pi/6) sin(6*pi/6) 0; cos(8*pi/6) sin(8*pi/6) 0; cos(10*pi/6) sin(10*pi/6) 0;
        1 0 0; cos(2*pi/6) 0 sin(2*pi/6); cos(4*pi/6) 0 sin(4*pi/6); cos(6*pi/6) 0 sin(6*pi/6); cos(8*pi/6) 0 sin(8*pi/6); cos(10*pi/6) 0 sin(10*pi/6);
        cos(10*pi/6) 0 sin(10*pi/6); 0 sqrt(1/2) -sqrt(1/2); cos(4*pi/6) sin(4*pi/6) 0; cos(4*pi/6) 0 sin(4*pi/6); 0 -sqrt(1/2) sqrt(1/2); cos(10*pi/6) sin(10*pi/6) 0;
        cos(8*pi/6) 0 sin(8*pi/6); 0 sqrt(1/2) -sqrt(1/2); cos(2*pi/6) sin(2*pi/6) 0; cos(2*pi/6) 0 sin(2*pi/6); 0 -sqrt(1/2) sqrt(1/2); cos(8*pi/6) sin(8*pi/6) 0]'
    conf = doubleRings(p0, Weave)
    initialConfiguration = toArray(conf, Weave) #+ (randn(Float64, length(toArray(p0, Weave))) .- 0.5)*0.05
    #q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-12)
    plotWeaving(initialConfiguration, Weave; colorscheme=[:yellow,:cyan,:purple,:yellow,:cyan,:purple])
    q = computeOptimalWeaving(q, Weave; tol=1e-2)
    normalspace=evaluate.(differentiate(Weave.constraints, Weave.coordinateVariables), Weave.coordinateVariables=>q)
    display("Tangent direction")
    display(toMatrix(collect(nullspace(normalspace)[1:end,1]), Weave))
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:yellow,:cyan,:purple])
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

function test_sphere_dodecahedron()
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
    plotWeaving(q, Weave; colorscheme=[:yellow,  :cyan,  :purple2, :green2, :blue2, :magenta2])
    toSphere(q, Weave)
end

function test_sphere_truncatedcube()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true], [(1, 40), (2, 39), (3, 16), (4, 15), (5, 43), (6, 42), (7, 28), (8, 27), (9, 38), (10, 37), (11, 24), (12, 23), (13, 45), (14, 44), (17, 36), (18, 35), (19, 32), (20, 31), (21, 47), (22, 46), (25, 34), (26, 33), (29, 41), (30, 48)], 
                                            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,1), (9,10),(10,11),(11,12),(12,13),(13,14),(14,15),(15,16),(16,9), (17,18),(18,19),(19,20),(20,21),(21,22),(22,23),(23,24),(24,17), (25,26),(26,27),(27,28),(28,29),(29,30),(30,31),(31,32),(32,25), (33,34),(34,35),(35,36),(36,37),(37,38),(38,39),(39,40),(40,33), (41,42),(42,43),(43,44),(44,45),(45,46),(46,47),(47,48),(48,41)], 
                                            []#=[(1,2,3,4,5,6,7,8,9,10), (11,12,13,14,15,16,17,18,19,20), (21,22,23,24,25,26,27,28,29,30), (31,32,33,34,35,36,37,38,39,40), (41,42,43,44,45,46,47,48,49,50), (51,52,53,54,55,56,57,58,59,60)]=#; samples = [[sqrt(2)-1,1,1]], offset=0.075)
    ξ = sqrt(2)-1
    p0 = [ξ 1 1; -ξ 1 1; -1 ξ 1; -1 -ξ 1; -ξ -1 1; ξ -1 1; 1 -ξ 1; 1 ξ 1;
        -1 1 ξ; -1 1 -ξ; -1 ξ -1; -1 -ξ -1; -1 -1 -ξ; -1 -1 ξ; -1 -ξ 1; -1 ξ 1;
        -ξ 1 -1; ξ 1 -1; 1 ξ -1; 1 -ξ -1; ξ -1 -1; -ξ -1 -1; -1 -ξ -1; -1 ξ -1;
        1 1 -ξ; 1 1 ξ; 1 ξ 1; 1 -ξ 1; 1 -1 ξ; 1 -1 -ξ; 1 -ξ -1; 1 ξ -1;
        1 1 ξ; 1 1 -ξ; ξ 1 -1; -ξ 1 -1; -1 1 -ξ; -1 1 ξ; -ξ 1 1; ξ 1 1;
        1 -1 ξ; ξ -1 1; -ξ -1 1; -1 -1 ξ; -1 -1 -ξ; -ξ -1 -1; ξ -1 -1; 1 -1 -ξ;
    ]' ./ sqrt(ξ^2+2)

    #=bars = []
    for i in 1:size(p0)[2]
        for j in i+1:size(p0)[2]
            if isapprox(norm(p0[:,i]-p0[:,j]),0)
                push!(bars,(i,j))
            end
        end
    end
    println(bars)=#

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-10)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:green,:teal,:red])
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:green,:teal,:red])
    toSphere(q, Weave)
end


function test_sphere_trunctetra()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false, false,true,false,true,false,true, false,true,false,true,false,true, false,true,false,true,false,true], [(1,13),(2,8),(3,9),(4,24),(5,19),(6,18),(7,14),(10,23),(11,22),(12,15),(16,21),(17,20)], 
                                            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1), (7,8),(8,9),(9,10),(10,11),(11,12),(12,7), (13,14),(14,15),(15,16),(16,17),(17,18),(18,13), (19,20),(20,21),(21,22),(22,23),(23,24),(24,19)], 
                                            []#=[(1,2,3,4,5,6,7,8,9,10), (11,12,13,14,15,16,17,18,19,20), (21,22,23,24,25,26,27,28,29,30), (31,32,33,34,35,36,37,38,39,40), (41,42,43,44,45,46,47,48,49,50), (51,52,53,54,55,56,57,58,59,60)]=#; offset=0.075)
    p0 = [1 3 1; 1 1 3; -1 -1 3; -3 -1 1; -3 1 -1; -1 3 -1; 
          3 1 1; 1 1 3; -1 -1 3; -1 -3 1; 1 -3 -1; 3 -1 -1; 
          1 3 1; 3 1 1; 3 -1 -1; 1 -1 -3; -1 1 -3; -1 3 -1;
          -3 1 -1; -1 1 -3; 1 -1 -3; 1 -3 -1; -1 -3 1; -3 -1 1;
    ]' ./ sqrt(11)

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple, :yellow,:cyan,:purple])
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple, :yellow,:cyan,:purple])
    toSphere(q, Weave)
end


function test_sphere_antiprism_square()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false], [(1,6), (2,13), (3,8), (4,15), (5,10), (7,12), (9,14), (11,16)], 
                                            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(10,11),(11,12),(12,13),(13,14),(14,15),(15,16),(16,1)], 
                                            []#=[(1,2,3,4,5,6,7,8,9,10), (11,12,13,14,15,16,17,18,19,20), (21,22,23,24,25,26,27,28,29,30), (31,32,33,34,35,36,37,38,39,40), (41,42,43,44,45,46,47,48,49,50), (51,52,53,54,55,56,57,58,59,60)]=#; samples=[[1,0,1]], offset=0.1)
    p0 = [1 0 1; 0 1 1; -1/sqrt(2) sqrt(1/2) -1; -1/sqrt(2) -sqrt(1/2) -1; 0 -1 1; 1 0 1;  1/sqrt(2) sqrt(1/2) -1; -1/sqrt(2) sqrt(1/2) -1;
         -1 0 1; 0 -1 1; 1/sqrt(2) -sqrt(1/2) -1; 1/sqrt(2) sqrt(1/2) -1; 0 1 1; -1 0 1; -1/sqrt(2) -sqrt(1/2) -1; 1/sqrt(2) -sqrt(1/2) -1;
    ]' ./ sqrt(2)

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple])
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple])
    toSphere(q, Weave)
end

function test_sphere_antiprism_penta()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false], [(1, 14), (2, 9), (3, 16), (4, 11), (5, 18), (6, 13), (7, 20), (8, 15), (10, 17), (12, 19)], 
                                            [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(10,11),(11,12),(12,13),(13,14),(14,15),(15,16),(16,17),(17,18),(18,19),(19,20),(20,1)], 
                                            []#=[(1,2,3,4,5,6,7,8,9,10), (11,12,13,14,15,16,17,18,19,20), (21,22,23,24,25,26,27,28,29,30), (31,32,33,34,35,36,37,38,39,40), (41,42,43,44,45,46,47,48,49,50), (51,52,53,54,55,56,57,58,59,60)]=#; samples=[], offset=0.075)
    p0 = [cos(0) sin(0) 1; cos(2*pi/5) sin(2*pi/5) 1; cos(3*pi/5) sin(3*pi/5) -1; cos(5*pi/5) sin(5*pi/5) -1; cos(6*pi/5) sin(6*pi/5) 1; cos(8*pi/5) sin(8*pi/5) 1; cos(9*pi/5) sin(9*pi/5) -1; cos(pi/5) sin(pi/5) -1;
            cos(2*pi/5) sin(2*pi/5) 1; cos(4*pi/5) sin(4*pi/5) 1; cos(5*pi/5) sin(5*pi/5) -1; cos(7*pi/5) sin(7*pi/5) -1; cos(8*pi/5) sin(8*pi/5) 1; cos(0*pi/5) sin(0*pi/5) 1; cos(1*pi/5) sin(1*pi/5) -1; cos(3*pi/5) sin(3*pi/5) -1; 
            cos(4*pi/5) sin(4*pi/5) 1; cos(6*pi/5) sin(6*pi/5) 1; cos(7*pi/5) sin(7*pi/5) -1; cos(9*pi/5) sin(9*pi/5) -1;
    ]' ./ sqrt(2)
    #=
    bars = []
    for i in 1:size(p0)[2]
        for j in i+1:size(p0)[2]
            if isapprox(norm(p0[:,i]-p0[:,j]),0)
                push!(bars,(i,j))
            end
        end
    end
    println(bars)=#

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple])
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple])
    toSphere(q, Weave)
end



function test_flattorus_plain()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true, false,true,false, true,false,true, false,true,false, true,false,true, false,true,false], [(1,10),(2,13),(3,16), (4,11),(5,14),(6,17), (7,12),(8,15),(9,18)], 
                                            [(1,2),(2,3), (4,5),(5,6), (7,8),(8,9), (10,11),(11,12), (13,14),(14,15), (16,17), (17,18)], 
                                            []; offset=0.1, manifold = "flattorus", samples = [[0.,0.,0],[]], PBC=[(3,1,[1.,0,0]), (6,4,[1.,0,0]), (9,7,[1.,0,0]), (12,10,[0.,1.,0]), (15,13,[0.,1.,0]), (18,16,[0.,1.,0])])
    p0 = [0 0 0; 0.5 0 0; 1 0 0; 0 0.5 0; 0.5 0.5 0; 1 0.5 0; 0 1 0; 0.5 1 0; 1 1 0; 
        0 0 0; 0 0.5 0; 0 1 0; 0.5 0 0; 0.5 0.5 0; 0.5 1 0; 1 0 0; 1 0.5 0; 1 1 0; 
    ]'

    initialConfiguration = toArray(p0,Weave)

    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave)

    q = computeOptimalWeaving(q, Weave; maxseconds=50)
    plotWeaving(q, Weave)
    toFlatTorus(q,Weave)

end

function test_torus_a()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true,false,true,false,true,false, false,true,false,true,false,true,false,true], [(1,9),(3,11),(5,13),(7,15)], [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,1), (9,10),(10,11),(11,12),(12,13),(13,14),(14,15),(15,16),(16,9)], []; samples=[], manifold="torus", offset=0.075)
    p0 = [-1.5 0 0; -1/sqrt(2) -1/sqrt(2) 0.5; 0 -0.5 0; 1/sqrt(2) -1/sqrt(2) -0.5; 1.5 0 0; 1/sqrt(2) 1/sqrt(2) 0.5; 0 0.5 0; -1/sqrt(2) 1/sqrt(2) -0.5;
    -1.5 0 0; -1/sqrt(2) -1/sqrt(2) -0.5; 0 -0.5 0; 1/sqrt(2) -1/sqrt(2) 0.5; 1.5 0 0; 1/sqrt(2) 1/sqrt(2) -0.5; 0 0.5 0; -1/sqrt(2) 1/sqrt(2) 0.5;
    ]'
    initialConfiguration = toArray(p0, Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave)
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave)
end



function test_torus_trefoil()
    Weave = WeavingsOnManifolds.WeavingOnManifold([true, true, true, true, true, true, true, true, true, true, true, true], [(1,7), (2,8), (3,9), (4,10), (5,11), (6,12)], [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(10,11),(11,12),(12,1)], []; torusknot=true, manifold="torus", samples=[[1,0,0.5]], offset=0.)
    p0 = [1 0 0.5; 1.5*cos(pi/3) 1.5*sin(pi/3) 0; cos(2*pi/3) sin(2*pi/3) -0.5; 0.5*cos(3*pi/3) 0.5*sin(3*pi/3) 0; cos(4*pi/3) sin(4*pi/3) 0.5; 1.5*cos(5*pi/3) 1.5*sin(5*pi/3) 0;
    1 0 -0.5; 0.5*cos(pi/3) 0.5*sin(pi/3) 0; cos(2*pi/3) sin(2*pi/3) 0.5; 1.5*cos(3*pi/3) 1.5*sin(3*pi/3) 0; cos(4*pi/3) sin(4*pi/3) -0.5; 0.5*cos(5*pi/3) 0.5*sin(5*pi/3) 0; 
    ]'
    initialConfiguration = toArray(p0, Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    toTorus(q,Weave)
    plotWeaving(q, Weave)
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave)
    toTorus(q,Weave)
end

function test_hyperboloid()
    c = 1.0025
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,true,true,true, true,false,true, true,false,true, true,false,true, true,false,true], [(1,6), (2,9), (3,12), (4,15)], [(1,2),(2,3),(3,4),(4,1), (5,6),(6,7), (8,9),(9,10), (11,12),(12,13), (14,15),(15,16)], []; 
        pinned_vertices = Dict(5=>c.*[1, -1, 1], 7=>c.*[1, 1, -1], 8=>c.*[1, 1, 1], 10=>c.*[-1, 1, -1], 11=>c.*[-1, 1, 1], 13=>c.*[-1, -1, -1], 14=>c.*[-1, -1, 1], 16=>c.*[1, -1, -1]),
        manifold="hyperboloid", samples=[], offset=0.075)
    display(Weave.constraints)
    p0 = [1 0 0; 0 1 0; -1 0 0; 0 -1 0; 
            sqrt(2) -sqrt(2) 1; 1 0 0; sqrt(2) sqrt(2) -1; sqrt(2) sqrt(2) 1; 0 1 0; -sqrt(2) sqrt(2) -1;
            -sqrt(2) sqrt(2) 1; -1 0 0; -sqrt(2) -sqrt(2) -1; -sqrt(2) -sqrt(2) 1; 0 -1 0; sqrt(2) -sqrt(2) -1;
    ]'
    display(toMatrix(toArray(p0, Weave), Weave))
    initialConfiguration = toArray(p0, Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave)
    q = computeOptimalWeaving(q, Weave)
    display(toMatrix(q, Weave))
    toHyperboloid(q,Weave)
end

function test_hyperboloid2()
    c = 1.0025
    Weave = WeavingsOnManifolds.WeavingOnManifold([true,false,true, false,true,false], [(2,5)], [(1,2),(2,3), (4,5),(5,6)], []; 
        pinned_vertices = Dict(1=>c.*[1, 0, sqrt(2)], 3=>c.*[-1, 0, sqrt(2)], 4=>c.*[0, 1, sqrt(2)], 6=>c.*[0, -1, sqrt(2)]),
        manifold="hyperboloid", samples=[], offset=0.075)
    display(Weave.constraints)
    p0 = [1 0 sqrt(2); 0 0 1; -1 0 sqrt(2); 0 1 sqrt(2); 0 0 1; 0 -1 sqrt(2);
    ]'
    display(toMatrix(toArray(p0, Weave), Weave))
    initialConfiguration = toArray(p0, Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave)
    q = computeOptimalWeaving(q, Weave)
    display(toMatrix(q, Weave))
    toHyperboloid(q,Weave)
end

function test_cylinder_kagome()
    Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true,false,false,true,false,true,false, true,false,true,false,true,true,false,true,false,true, false,true,false,true,true,false,true,false], 
    [(1,20),(2,21), (3,13), (4,26), (5,16), (6,15), (7,23), (8,18), (9,28),(10,11), (12,25), (14,22), (17,27), (19,24)], [(1,2),(2,3),(3,4),(4,5), (6,7),(7,8),(8,9),(9,10), (11,12),(12,13),(13,14),(14,15), (16,17),(17,18),(18,19),(19,20), (21,22),(22,23),(23,24),(24,21), (25,26),(26,27),(27,28),(28,25)], []; 
        manifold="cylinder", samples=[], offset=0.075, PBC=[(1,10,[0,0,1]), (6,5,[0,0,1]), (15,16,[0,0,1]), (20,11,[0,0,1])],
        pinned_vertices = Dict(6=>[[3],[0.33333]], 15=>[[3],[0.66666]], 20=>[[3],[1]]))
    p0 = [1 0 0; cos(pi/4) sin(pi/4) 0.25; 0 1 0.5; cos(3*pi/4) sin(3*pi/4) 0.75; -1 0 1; -1 0 0; cos(5*pi/4) sin(5*pi/4) 0.25; 0 -1 0.5; cos(7*pi/4) sin(7*pi/4) 0.75; 1 0 1; 
        1 0 1; cos(pi/4) sin(pi/4) 0.75; 0 1 0.5; cos(3*pi/4) sin(3*pi/4) 0.25; -1 0 0; -1 0 1; cos(5*pi/4) sin(5*pi/4) 0.75; 0 -1 0.5; cos(7*pi/4) sin(7*pi/4) 0.25; 1 0 0; 
        cos(pi/4) sin(pi/4) 0.25; cos(3*pi/4) sin(3*pi/4) 0.25; cos(5*pi/4) sin(5*pi/4) 0.25; cos(7*pi/4) sin(7*pi/4) 0.25;
        cos(pi/4) sin(pi/4) 0.75; cos(3*pi/4) sin(3*pi/4) 0.75; cos(5*pi/4) sin(5*pi/4) 0.75; cos(7*pi/4) sin(7*pi/4) 0.75;
    ]'
    display(toMatrix(toArray(p0, Weave), Weave))
    initialConfiguration = toArray(p0, Weave)

    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    plotWeaving(q, Weave)
    println(toMatrix(q,Weave))

    q = computeOptimalWeaving(q, Weave)
    display(toMatrix(q, Weave))
    plotWeaving(q, Weave)
    toCylinder(q,Weave)
end

function test_cylinder_kagome2()
    Weave = WeavingsOnManifolds.WeavingOnManifold([false,true,false,true,false,false,true,false,true,false, true,false,true,false,true,true,false,true,false,true, false,true,false,true,false,false,true,false,true,false, true,false,true,false,true,true,false,true,false,true, false,true,false,true,false,true,false,true, true,false,true,false,true,false,true,false], 
    [(1, 40), (2, 41), (3, 13), (4, 50), (5, 16), (6, 15), (7, 43), (8, 18), (9, 52), (10, 31), (11,30), (12, 49), (14, 42), (17, 51), (19, 44), (20, 21), (22, 45), (23, 33), (24, 54), (25, 36), (26, 35), (27, 47), (28, 38), (29, 56), (32, 53), (34, 46), (37, 55), (39, 48)], 
    [(1,2),(2,3),(3,4),(4,5), (6,7),(7,8),(8,9),(9,10), (11,12),(12,13),(13,14),(14,15), (16,17),(17,18),(18,19),(19,20), (21,22),(22,23),(23,24),(24,25), (26,27),(27,28),(28,29),(29,30), (31,32),(32,33),(33,34),(34,35), (36,37),(37,38),(38,39),(39,40), (41,42),(42,43),(43,44),(44,45),(45,46),(46,47),(47,48),(48,41), (49,50),(50,51),(51,52),(52,53),(53,54),(54,55),(55,56),(56,49)], []; 
        manifold="cylinder", samples=[], offset=0.065, PBC=[(6, 5, [0, 0, 1]), (15, 16, [0, 0, 1]), (26, 25, [0, 0, 1]), (35, 36, [0, 0, 1]), (1, 30, [0,0,1]), (40, 11,[0,0,1]), (21,10, [0,0,1]), (20, 31, [0,0,1])],
        pinned_vertices = Dict())
    p0 = [1 0 0; cos(pi/8) sin(pi/8) 0.25; cos(2*pi/8) sin(2*pi/8) 0.5; cos(3*pi/8) sin(3*pi/8) 0.75; cos(4*pi/8) sin(4*pi/8) 1; cos(4*pi/8) sin(4*pi/8) 0; cos(5*pi/8) sin(5*pi/8) 0.25; cos(6*pi/8) sin(6*pi/8) 0.5; cos(7*pi/8) sin(7*pi/8) 0.75; cos(8*pi/8) sin(8*pi/8) 1; 
        1 0 1; cos(pi/8) sin(pi/8) 0.75; cos(2*pi/8) sin(2*pi/8) 0.5; cos(3*pi/8) sin(3*pi/8) 0.25; cos(4*pi/8) sin(4*pi/8) 0; cos(4*pi/8) sin(4*pi/8) 1; cos(5*pi/8) sin(5*pi/8) 0.75; cos(6*pi/8) sin(6*pi/8) 0.5; cos(7*pi/8) sin(7*pi/8) 0.25; cos(8*pi/8) sin(8*pi/8) 0; 

        cos(pi) sin(pi) 0; cos(pi+pi/8) sin(pi+pi/8) 0.25; cos(pi+2*pi/8) sin(pi+2*pi/8) 0.5; cos(pi+3*pi/8) sin(pi+3*pi/8) 0.75; cos(pi+4*pi/8) sin(pi+4*pi/8) 1; cos(pi+4*pi/8) sin(pi+4*pi/8) 0; cos(pi+5*pi/8) sin(pi+5*pi/8) 0.25; cos(pi+6*pi/8) sin(pi+6*pi/8) 0.5; cos(pi+7*pi/8) sin(pi+7*pi/8) 0.75; cos(pi+8*pi/8) sin(pi+8*pi/8) 1; 
        cos(pi) sin(pi) 1; cos(pi+pi/8) sin(pi+pi/8) 0.75; cos(pi+2*pi/8) sin(pi+2*pi/8) 0.5; cos(pi+3*pi/8) sin(pi+3*pi/8) 0.25; cos(pi+4*pi/8) sin(pi+4*pi/8) 0; cos(pi+4*pi/8) sin(pi+4*pi/8) 1; cos(pi+5*pi/8) sin(pi+5*pi/8) 0.75; cos(pi+6*pi/8) sin(pi+6*pi/8) 0.5; cos(pi+7*pi/8) sin(pi+7*pi/8) 0.25; cos(pi+8*pi/8) sin(pi+8*pi/8) 0; 

        cos(pi/8) sin(pi/8) 0.25; cos(3*pi/8) sin(3*pi/8) 0.25; cos(5*pi/8) sin(5*pi/8) 0.25; cos(7*pi/8) sin(7*pi/8) 0.25; cos(pi+pi/8) sin(pi+pi/8) 0.25; cos(pi+3*pi/8) sin(pi+3*pi/8) 0.25; cos(pi+5*pi/8) sin(pi+5*pi/8) 0.25; cos(pi+7*pi/8) sin(pi+7*pi/8) 0.25;
        cos(pi/8) sin(pi/8) 0.75; cos(3*pi/8) sin(3*pi/8) 0.75; cos(5*pi/8) sin(5*pi/8) 0.75; cos(7*pi/8) sin(7*pi/8) 0.75; cos(pi+pi/8) sin(pi+pi/8) 0.75; cos(pi+3*pi/8) sin(pi+3*pi/8) 0.75; cos(pi+5*pi/8) sin(pi+5*pi/8) 0.75; cos(pi+7*pi/8) sin(pi+7*pi/8) 0.75;
    ]'

    #=
    pbc = []
    for n in 1:size(p0)[2], m in n+1:size(p0)[2]
        if isapprox(norm(p0[:,n]-p0[:,m]-[0,0,1]), 0) && Weave.offsetList[n]==Weave.offsetList[m]
            push!(pbc, (m,n,[0,0,1]))
        elseif isapprox(norm(p0[:,n]-p0[:,m]+[0,0,1]), 0) && Weave.offsetList[n]==Weave.offsetList[m]
            push!(pbc, (n,m,[0,0,1]))
        end
    end
    println(pbc)=#

    println(toMatrix(toArray(p0, Weave), Weave))
    initialConfiguration = toArray(p0, Weave)
    display("text/plain",toMatrix(initialConfiguration,Weave))    
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-8)
    println(evaluate(Weave.constraints,Weave.coordinateVariables=>q))
    plotWeaving(q, Weave)
    display(toMatrix(q,Weave)')
    toCylinder(q,Weave)
    q = computeOptimalWeaving(q, Weave)
    display(toMatrix(q, Weave))
    plotWeaving(q, Weave)
    toCylinder(q,Weave)
end

function test_sphere_truncated_octahedron()
    offset = 0.065
    pointlist, edgelist = read_poly_file("trunc_oct.poly")
    newedgelist = [edgelist[1]]
    deleteat!(edgelist, 1)
    while ! isempty(edgelist)
        index = findfirst(ring -> any(edge -> (edge in ring)||((edge[2],edge[1]) in ring), vcat(newedgelist...)), edgelist)
        push!(newedgelist, popat!(edgelist, index))
    end
    cables = []
    for ring in newedgelist
        ring_sugg = [(ring[i][1]+mod(i,2)*length(pointlist), ring[i][2]+mod(i+1,2)*length(pointlist)) for i in 1:length(ring)]
        if !any(ringlist -> any(cable->ringlist==cable||ringlist==(cable[2],cable[1]), cables), ring_sugg)
            #println(filter(t->any(entry->entry in t, vcat([[r[1], r[2]] for r in ring_sugg]...)), cables), ring_sugg)
            append!(cables, ring_sugg)
        else
            ring_sugg = [(ring[i][1]+mod(i+1,2)*length(pointlist), ring[i][2]+mod(i,2)*length(pointlist)) for i in 1:length(ring)]
            if any(ringlist -> any(cable->ringlist==cable||ringlist==(cable[2],cable[1]), cables), ring_sugg)
                continue
            end    
            append!(cables, ring_sugg)
        end
        println("")
    end
    p0 = hcat(vcat(pointlist .* (1+offset), pointlist .* (1-offset))...)
    barlist = [(bar, bar+length(pointlist)) for bar in 1:length(pointlist)]
    offsetlist = vcat([true for _ in 1:length(pointlist)], [false for _ in 1:length(pointlist)])

    Weave = WeavingsOnManifolds.WeavingOnManifold(offsetlist, barlist, cables, []; samples = [pointlist[1]], offset=offset)

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-10)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:green,:teal,:red,:yellow,:cyan,:purple,:green,:teal,:red])
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:green,:teal,:red,:yellow,:cyan,:purple,:green,:teal,:red])
    toSphere(q, Weave)
end


function test_sphere_truncated_icosihedron()
    offset = 0.06
    pointlist, edgelist = read_poly_file("trunc_icos.poly")
    newedgelist = [edgelist[1]]
    deleteat!(edgelist, 1)
    while ! isempty(edgelist)
        index = findfirst(ring -> any(edge -> (edge in ring)||((edge[2],edge[1]) in ring), vcat(newedgelist...)), edgelist)
        push!(newedgelist, popat!(edgelist, index))
    end
    cables = []
    for ring in newedgelist
        ring_sugg = [(ring[i][1]+mod(i,2)*length(pointlist), ring[i][2]+mod(i+1,2)*length(pointlist)) for i in 1:length(ring)]
        if !any(ringlist -> any(cable->ringlist==cable||ringlist==(cable[2],cable[1]), cables), ring_sugg)
            #println(filter(t->any(entry->entry in t, vcat([[r[1], r[2]] for r in ring_sugg]...)), cables), ring_sugg)
            append!(cables, ring_sugg)
        else
            ring_sugg = [(ring[i][1]+mod(i+1,2)*length(pointlist), ring[i][2]+mod(i,2)*length(pointlist)) for i in 1:length(ring)]
            if any(ringlist -> any(cable->ringlist==cable||ringlist==(cable[2],cable[1]), cables), ring_sugg)
                continue
            end    
            append!(cables, ring_sugg)
        end
        println("")
    end
    p0 = hcat(vcat(pointlist .* (1+offset), pointlist .* (1-offset))...)
    barlist = [(bar, bar+length(pointlist)) for bar in 1:length(pointlist)]
    offsetlist = vcat([true for _ in 1:length(pointlist)], [false for _ in 1:length(pointlist)])

    Weave = WeavingsOnManifolds.WeavingOnManifold(offsetlist, barlist, cables, []; samples = [pointlist[1]], offset=offset)

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-10)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:green,:teal,:red,:yellow,:cyan,:purple,:green,:teal,:red,:yellow,:cyan,:purple,:green,:teal,:red,:yellow,:cyan,:purple,:green,:teal,:red])
    q = computeOptimalWeaving(q, Weave)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:green,:teal,:red,:yellow,:cyan,:purple,:green,:teal,:red,:yellow,:cyan,:purple,:green,:teal,:red,:yellow,:cyan,:purple,:green,:teal,:red])
    toSphere(q, Weave)

end


function test_sphere_truncated_dodecahedron()
    offset = 0.06
    pointlist, edgelist = read_poly_file("trunc_dodec.poly")
    newedgelist = [edgelist[1]]
    deleteat!(edgelist, 1)
    while ! isempty(edgelist)
        index = findfirst(ring -> any(edge -> (edge in ring)||((edge[2],edge[1]) in ring), vcat(newedgelist...)), edgelist)
        push!(newedgelist, popat!(edgelist, index))
    end
    cables = []
    for ring in newedgelist
        ring_sugg = [(ring[i][1]+mod(i,2)*length(pointlist), ring[i][2]+mod(i+1,2)*length(pointlist)) for i in 1:length(ring)]
        if !any(ringlist -> any(cable->ringlist==cable||ringlist==(cable[2],cable[1]), cables), ring_sugg)
            #println(filter(t->any(entry->entry in t, vcat([[r[1], r[2]] for r in ring_sugg]...)), cables), ring_sugg)
            append!(cables, ring_sugg)
        else
            ring_sugg = [(ring[i][1]+mod(i+1,2)*length(pointlist), ring[i][2]+mod(i,2)*length(pointlist)) for i in 1:length(ring)]
            if any(ringlist -> any(cable->ringlist==cable||ringlist==(cable[2],cable[1]), cables), ring_sugg)
                continue
            end
            append!(cables, ring_sugg)
        end
    end
    p0 = hcat(vcat(pointlist .* (1+offset), pointlist .* (1-offset))...)
    barlist = [(bar, bar+length(pointlist)) for bar in 1:length(pointlist)]
    offsetlist = vcat([true for _ in 1:length(pointlist)], [false for _ in 1:length(pointlist)])

    Weave = WeavingsOnManifolds.WeavingOnManifold(offsetlist, barlist, cables, []; samples = [pointlist[1] .* (1+offset)], offset=offset)

    initialConfiguration = toArray(p0,Weave)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints; tol = 1e-10)
    toSphere(q, Weave)
    plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:green,:teal,:red, :yellow,:cyan,:purple,:green,:teal,:red])
    for _ in 1:25
        q = computeOptimalWeaving(q, Weave)
        plotWeaving(q, Weave; colorscheme=[:yellow,:cyan,:purple,:green,:teal,:red, :yellow,:cyan,:purple,:green,:teal,:red])
        toSphere(q, Weave)
    end
end


function newtonCorrect(q::Vector{Float64}, variables, equations; tol = 1e-8)
    global qnew, damping = q, 0.25
	jac = hcat([differentiate(eq, variables) for eq in equations]...)
    while( norm(evaluate.(equations, variables=>q)) > tol )
        println(norm(evaluate.(equations, variables=>q)))
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

function read_poly_file(str::String)
    pointlist, edgelist = [], []
    open(str, "r") do file
        s = readline(file)
        while ! eof(file) 
              s = readline(file)
              if occursin("POLYS",s)
                break
              end
              new_s = split(s,": ")[2]
              point = parse.(Float64, split(new_s, " "))
              push!(pointlist, point)
        end

        while ! eof(file) 
            s = readline(file)
            if occursin("END",s)
              break
            end
            new_s = split(s,": ")[2]
            edges = parse.(Int64, split(new_s, " "))
            helper_edges=[]
            for i in 1:length(edges)-1
                push!(helper_edges, (edges[i],edges[i+1]))
            end
            push!(edgelist, helper_edges)
        end
    end
    return pointlist, edgelist
end

function write_txt_file(energyArray; filename)
    open("$(filename)EnergyFunctionals.txt", "w") do file
        for ar in energyArray 
            write(file, "$(ar[1]) $(ar[2]) $(ar[3]) $(ar[4])\n")
        end
    end
end


function computeOptimalWeaving(initialConfiguration::Vector{Float64}, Weave::WeavingOnManifold; tol = 1e-3, maxseconds=600, bending = bending)
    q = newtonCorrect(initialConfiguration, Weave.coordinateVariables, Weave.constraints)
    Q = x->energyFunction(x, Weave; bending = bending)
    G = HomotopyOpt.ConstraintVariety(Weave.coordinateVariables, Weave.constraints, length(Weave.coordinateVariables), length(Weave.coordinateVariables)-length(Weave.constraints))
    #@time resultmin = findminima(q, tol, G, Q; whichstep="gaussnewtonstep", stepdirection = "gradientdescent", maxseconds=maxseconds)
    #@time resultmin = HomotopyOpt.findminima(q, tol, G, Q; whichstep="HomotopyContinuation", stepdirection = "gradientdescent", maxseconds=maxseconds)
    #@time resultmin = HomotopyOpt.findminima(q, tol, G, Q; whichstep="EDStep", stepdirection = "gradientdescent", maxseconds=maxseconds)
    @time resultmin = HomotopyOpt.findminima(q, tol, G, Q; whichstep="Algorithm 0", homotopyMethod = "HomotopyContinuation", maxseconds=maxseconds, initialstepsize=0.01)
    return(resultmin.computedpoints[end])
end

function plot_energy(EnergyArray, Weave; pltlabel="", stepsize=0.01, endpoint=0.5, bending = false)
    plt = Plots.plot(1:length(EnergyArray), [en[1]/EnergyArray[1][1] for en in EnergyArray], xticks = (1:50:201, string.([1,0.75,0.5,0.75,1])); size=(900,700), xtickfontsize=15,ytickfontsize=15, legendfontsize=15, linewidth=6, color=:steelblue, label = "Total Energy", grid=false, legend_column=bending ? 2 : 1)
    Plots.plot!(plt, 1:length(EnergyArray), [en[2]/EnergyArray[1][1] for en in EnergyArray], xticks = (1:50:201, string.([1,0.75,0.5,0.75,1])); linewidth=6, xtickfontsize=15,ytickfontsize=15, legendfontsize=15, color=:forestgreen, label = "Cables")
    Plots.plot!(plt, 1:length(EnergyArray), [en[3]/EnergyArray[1][1] for en in EnergyArray], xticks = (1:50:201, string.([1,0.75,0.5,0.75,1])); linewidth=6, xtickfontsize=15,ytickfontsize=15, legendfontsize=15, color=:firebrick2, label = "Struts")
    if bending
        Plots.plot!(plt, 1:length(EnergyArray), [en[4]/EnergyArray[1][1] for en in EnergyArray], xticks = (1:50:201, string.([1,0.75,0.5,0.75,1])); linewidth=6, xtickfontsize=15,ytickfontsize=15, legendfontsize=15, color=:grey, label = "Bending")
    end
    Plots.ylims!(plt,(minimum(vcat([en[2]/EnergyArray[1][1] for en in EnergyArray],[en[3]/EnergyArray[1][1] for en in EnergyArray],[en[4]/EnergyArray[1][1] for en in EnergyArray]))-0.05,maximum([en[1]/EnergyArray[1][1] for en in EnergyArray])+0.25))
    Plots.savefig(plt, "$(pltlabel)EnergyPlot.png")
end

function plot_hystereces(; prefix="src/")
    namespace = ["$(prefix)Tetrahedron3EnergyFunctionals.txt", "$(prefix)Octahedron3EnergyFunctionals.txt", "$(prefix)Dodecahedron3EnergyFunctionals.txt", "$(prefix)RhombiOctahedron3EnergyFunctionals.txt", "$(prefix)TruncTetra3EnergyFunctionals.txt", "$(prefix)TruncCube3EnergyFunctionals.txt", "$(prefix)TruncOct3EnergyFunctionals.txt", "$(prefix)TruncIcos3EnergyFunctionals.txt", "$(prefix)TruncDod3EnergyFunctionals.txt"]
    colorspace = [:yellow3, :blue2, :green2, :purple2, :red2, :cyan2, :magenta2, :orange2, :teal2]
    colorspace = Plots.distinguishable_colors(10, [Plots.RGB(1,1,1), Plots.RGB(0,0,0)], dropseed=false;)[2:end]
    labelspace= [L"\,[1/2]_{tet}", L"\,[1/2]_{cub}", L"\,[1/2]_{dod}", L"\,Rhombi", L"\,[2/2]_{tet}", L"\,[2/2]_{cub}", L"\,[2/2]_{oct}", L"\,[2/2]_{icos}", L"\,[2/2]_{dodec}"]
    EnergyArray = []
    for filename in namespace
        try
            helper = []
            display(filename)
            open("$(filename)", "r") do f
                while !eof(f)
                    s = readline(f)
                    if s==""
                        continue
                    end
                    s = split(s, " ")
                    vals = parse.( Float64, s)[1:4]
                    push!(helper, vals)
                end
            end
            push!(EnergyArray, helper)
        catch e
            display(e)
            nothing
        end
    end
    for en in EnergyArray
        display(en)
    end
    EnergyArray[3][101:201] = EnergyArray[3][101:-1:1]
    EnergyArray[6][1:101] = EnergyArray[6][201:-1:101]
    EnergyArray[7][1:101] = EnergyArray[7][201:-1:101]
    EnergyArray[8] = vcat(EnergyArray[8][1:101], EnergyArray[8][100:-1:1])
    EnergyArray[9] = vcat(EnergyArray[9][1:101], EnergyArray[9][100:-1:1])

    plt1 = Plots.plot(1:length(EnergyArray[1]), [(en[1]-en[4])/(EnergyArray[1][1][1]-EnergyArray[1][1][4]) for en in EnergyArray[1]], xticks = (1:50:201, string.([1,1.5,2,1.5,1])); size=(900,700), xtickfontsize=14,ytickfontsize=14, legendfontsize=18, titlefontsize=16, linewidth=5, title="Total Energy", color=colorspace[1], label = labelspace[1], grid=false, legend_column= 3, legend=:top)
    foreach(i->Plots.plot!(plt1, 1:length(EnergyArray[i]), [(en[1]-en[4])/(EnergyArray[i][1][1]-EnergyArray[i][1][4]) for en in EnergyArray[i]], xticks = (1:50:201, string.([1,1.5,2,1.5,1])); size=(900,700), xtickfontsize=14,ytickfontsize=14, legendfontsize=18, linewidth=5, color=colorspace[i], label = labelspace[i]), 2:length(EnergyArray))
    Plots.ylims!(plt1,(0.98725,1.0282))
    plt2 = Plots.plot(1:length(EnergyArray[1]), [en[2]/EnergyArray[1][1][2] for en in EnergyArray[1]], xticks = (1:50:201, string.([1,1.5,2,1.5,1])); size=(900,700), xtickfontsize=14,ytickfontsize=14, legendfontsize=16, titlefontsize=16, linewidth=5, title="Cable Energy", color=colorspace[1], label = labelspace[1], grid=false, legend=false)
    foreach(i->Plots.plot!(plt2, 1:length(EnergyArray[i]), [en[2]/EnergyArray[i][1][2] for en in EnergyArray[i]], xticks = (1:50:201, string.([1,1.5,2,1.5,1])); size=(900,700), xtickfontsize=14,ytickfontsize=14, legendfontsize=14, linewidth=5, color=colorspace[i], label = labelspace[i]), 2:length(EnergyArray))
    Plots.ylims!(plt2,(0.95, 1.13))
    plt = Plots.plot(plt1,plt2,layout=(1,2), size=(1400,600))
    Plots.savefig(plt, "HysterecisEnergyPlots_squashedsphere.png")
    display(plt)
end

display("start")
#test_spheroid_tetrahedron()
#test_spheroid_octahedron()
#test_spheroid_dodecahedron()
#test_spheroid_rhombioctahedron()
#test_spheroid_trunctetra()
#test_spheroid_truncatedcube()
#test_spheroid_truncatedoctahedron()
#test_spheroid_truncatedicosihedron()
#test_spheroid_truncateddodecahedron()
plot_hystereces()

#TODO how does the energy behave, when I squish the sphere => become an ellipse. 
#TODO how to restrict weavings on manifolds?
#TODO Tensegrities in the plane with saddle?
#TODO What if the cables are under compression instead of tension?W
#TODO What is the optimal manifold for the weaving => random 4-coordinated graph, alternating. Fixed volume.
#TODO Hyperauxeticity paper
end
