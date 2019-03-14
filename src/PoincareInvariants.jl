__precompile__()

module PoincareInvariants

    using ApproxFun
    using HDF5
    using ProgressMeter
    using GeometricIntegrators
    using GeometricIntegrators.Utils

    export PoincareInvariant1st, PoincareInvariant1stCanonical,
           PoincareInvariant2nd, PoincareInvariant2ndCanonical,
           PoincareInvariant2ndTrapezoidal,
           evaluate_poincare_invariant, write_to_hdf5

    include("matrix_utils.jl")
    include("poincare_invariant_1st_common.jl")
    include("poincare_invariant_1st.jl")
    include("poincare_invariant_1st_canonical.jl")
    include("poincare_invariant_2nd.jl")
    include("poincare_invariant_2nd_canonical.jl")
    include("poincare_invariant_2nd_trapezoidal.jl")
    include("poincare_invariant_2nd_common.jl")

end
