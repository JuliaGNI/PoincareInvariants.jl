# Reference

```@meta
DocTestSetup  = quote
    using PoincareInvariants
end
```

```@docs
PoincareInvariants
```

## Fundamentals

```@docs
compute!
FirstPoincareInvariant
FirstPoincareInvariant{T, D}(θ::θT, N::Integer, plan::P) where {T, D, θT, P}
FirstPoincareInvariant{T, D}(θ::θT, N::Integer, P::Type=DEFAULT_FIRST_PLAN) where {T, D, θT}
FirstPoincareInvariant{T, D, typeof(canonical_one_form)}(N::Integer, P=DEFAULT_FIRST_PLAN) where {T, D}
SecondPoincareInvariant
SecondPoincareInvariant{T, D}(ω::ωT, N, plan::P) where {T, D, ωT, P}
SecondPoincareInvariant{T, D}(ω, N, P::Type=DEFAULT_SECOND_PLAN) where {T, D}
SecondPoincareInvariant{T, D, CanonicalSymplecticMatrix{T}}(N, P::Type=DEFAULT_SECOND_PLAN) where {T, D}
```

## Interface

```@docs
AbstractPoincareInvariant
getpoints
getpointnum
getpointspec
getdim
getform
getplan
```

## Canonical Symplectic Forms

```@docs
CanonicalSymplecticMatrix
```
