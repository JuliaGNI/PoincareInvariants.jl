# Canonical Symplectic Structures

```@meta
CurrentModule = PoincareInvariants.CanonicalSymplecticStructures
DocTestSetup = quote
    using PoincareInvariants.CanonicalSymplecticStructures
end
```

This module currently contains the type [`CanonicalSymplecticTwoForm`](@ref), which satisifes the `AbstractMatrix` interface and has an optimised `LinearAlgebra.dot` method to compute the vector matrix vector product.

```@docs
CanonicalSymplecticTwoForm
```
