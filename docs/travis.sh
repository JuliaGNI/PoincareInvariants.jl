#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia --color=yes -e 'using Pkg; Pkg.add("Documenter"); Pkg.instantiate(); import PoincareInvariants; include(joinpath(dirname(pathof(PoincareInvariants)), "..", "docs", "make.jl"))';
fi
