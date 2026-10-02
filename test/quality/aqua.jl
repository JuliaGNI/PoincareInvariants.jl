using Aqua
using PoincareInvariants

# On Julia 1.11 the persistent-tasks wrapper's `Pkg.precompile` also builds
# `PoincareInvariantsMakieExt` after the package loads. That took 58 s and 97 s in two runs,
# more than Aqua's default `tmax` of 30 s. A real persistent task blocks forever, so a larger
# `tmax` hides none.
Aqua.test_all(PoincareInvariants; persistent_tasks = (; tmax = 300))
