module PoincareInvariants

    using HDF5
    using FFTW
    using OffsetArrays
    using ProgressMeter

    using GeometricIntegrators
    using GeometricIntegrators.Utils

    import DoubleFloats: Double64
    import LinearAlgebra: ⋅

    include("matrix_utils.jl")

    export PoincareInvariant1st, PoincareInvariant1stCanonical,
           PoincareInvariant2nd, PoincareInvariant2ndCanonical,
           PoincareInvariant2ndApproxFun,
           PoincareInvariant2ndApproxFunCanonical,
           PoincareInvariant2ndOPQ,
           PoincareInvariant2ndTrapezoidal,
           evaluate_poincare_invariant,
           evaluate_poincare_invariant_correction,
           write_to_hdf5

    abstract type AbstractPoincareInvariant2nd{DT} end
    abstract type AbstractPoincareInvariant2ndApproxFun{DT,ND,NC,NV} <: AbstractPoincareInvariant2nd{DT} end
    
    include("poincare_invariant_1st_common.jl")
    include("poincare_invariant_1st.jl")
    include("poincare_invariant_1st_canonical.jl")
    include("poincare_invariant_2nd_approxfun.jl")
    include("poincare_invariant_2nd_opq.jl")
    include("poincare_invariant_2nd_trapezoidal.jl")

    const PoincareInvariant2nd = PoincareInvariant2ndApproxFun
    const PoincareInvariant2ndCanonical = PoincareInvariant2ndApproxFunCanonical

end
