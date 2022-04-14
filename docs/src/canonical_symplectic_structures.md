# Canonical Symplectic Structures

```@meta
CurrentModule = PoincareInvariants.CanonicalSymplecticForms
DocTestSetup = quote
    using PoincareInvariants.CanonicalSymplecticForms
end
```

This module currently contains the type [`CanonicalSymplecticMatrix`](@ref), which satisifes the `AbstractMatrix` interface and has an optimised `LinearAlgebra.dot` method to compute the vector matrix vector product.

```@docs
CanonicalSymplecticMatrix
```
