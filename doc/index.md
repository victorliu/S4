# S4: Stanford Stratified Structure Solver

A program for computing electromagnetic fields in periodic, layered structures,
developed by Victor Liu (victorliu@alumni.stanford.edu),
graduate of the Fan group in the Stanford Electrical Engineering Department.

S4 consists of a core library with a C API that allows creating, specifying, and querying simulations.
It also contains Lua and Python frontend modules that can be imported into the corresponding base languages.
Finally, a separate standalone embedded Lua interpreter is provided for running simulations on massively parallel computing architectures supporting MPI or POSIX threads.

The Lua API (version 2) is [documented here](S4v2lua.md).
The base C API is [documented here](Capi.md).
