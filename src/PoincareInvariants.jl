module PoincareInvariants

    using HDF5
    using FFTW
    using OffsetArrays
    using ProgressMeter
    using Requires

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
    include("poincare_invariant_2nd_trapezoidal.jl")

    @require OrthogonalPolynomialsQuasi = "aa41a628-2c43-45df-899b-83ab96621781" begin
        @require FastTransforms = "057dd010-8810-581a-b7be-e3fc3b93f78c" begin
            include("poincare_invariant_2nd_opq.jl")
        end
    end

    const PoincareInvariant2nd = PoincareInvariant2ndApproxFun
    const PoincareInvariant2ndCanonical = PoincareInvariant2ndApproxFunCanonical

end
