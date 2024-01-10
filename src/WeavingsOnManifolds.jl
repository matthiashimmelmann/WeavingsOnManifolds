module WeavingsOnManifolds
import HomotopyContinuation: Expression, Variable, @var
import LinearAlgebra
import GLMakie

mutable struct WeavingOnManifold
    constraints
    energyFunction
    lagrangeSystem
    coordinateVariables
    lagrangeMultipliers

    function WeavingOnManifold(offsetList::Vector{Bool}, bars::Vector, cables::Vector, ambientVariables::Vector{Variable}, manifold::String, offset::Float64)
        (offset<0.5 && offset>0.05) || warn("We recommend an offset between 5% and 50%.")
        @var xs[1:3, 1:length(offsetList)]
        manifoldEquations, barEquations = [], []
        if isequal(manifold, "Sphere")||isequal(manifold, "sphere")
            manifoldEquations = [offsetList[i] ? sum( xs[:,i].^2 ) - (1+offset)^2 : sum( xs[:,i].^2 ) - (1-offset)^2 for i in 1:length(offsetList)]
        else
            throw(error("Only the Sphere is implemented right now!"))
        end
        barEquations = [sum( (xs[:,bar[1]] - xs[:,bar[2]]).^2 ) - (2*offset)^2 for bar in bars]
        energyFunction = sum([sum((xs[:,cable[1]] - xs[:,cable[2]]).^2) for cable in cables])
        @var λ[1:(length(manifoldEquations)+length(barEquations))]
        new(vcat(manifoldEquations, barEquations), energyFunction, energyFunction + λ'*vcat(manifoldEquations, barEquations), reduce(vcat, xs), λ)
    end
end

function computeOptimalWeaving(, Weaving::WeavingOnManifold)

end

end
