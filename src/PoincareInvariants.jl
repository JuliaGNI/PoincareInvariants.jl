module PoincareInvariants

# calculates phase points using param_func and outputs a vector of values for each
# dimension of phase space (as opposed to vector of SVectors)
#
# Type annotations added to force specialisation on function
#
# TODO: add @generated version with @nexprs
# TODO: Should generated version also use known point number?
function get_phase_points!(param_func::F, out, param_points, ::Val{N}) where {F, N}
    N::Int  # base Julia ntuple.jl does this, too
    
    for (i, pa_pnt) in enumerate(param_points)
        ph_pnt = param_func(pa_pnt)
        for j in 1:N
            out[j][i] = ph_pnt[j]
        end
    end

    out
end

end  # module
