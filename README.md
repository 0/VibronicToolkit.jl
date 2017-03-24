# VibronicToolkit.jl

Tools for vibronic Hamiltonians.
This package aims for correctness and readability rather than performance, so these tools are mainly useful for verifying other implementations.

Tested with Julia 0.5.


## Installation

1. `Pkg.clone("https://github.com/0/VibronicToolkit.jl.git")`


### Requirements

These should be pulled in automatically when installing this package.
To use it without installing it (e.g. from a local git checkout), you'll need to manually obtain the following dependencies:

* ArgParse (`Pkg.add("ArgParse")`)
* Distributions (`Pkg.add("Distributions")`)
* JSON (`Pkg.add("JSON")`)
* ProgressMeter (`Pkg.add("ProgressMeter")`)
* StatsBase (`Pkg.add("StatsBase")`)


## Units

* The input energy parameters define the unit of energy `E`.
* The input frequency parameters are in units of `E` (i.e. `hbar = 1`).
* Positions are dimensionless, so the other input parameters are also in units of `E`.
* Reciprocal temperatures are in units of reciprocal energy (i.e. `k_B = 1`).
* Heat capacities are in units of `E^2 / k_B T^2`.


## Examples

* `julia bin/analytical.jl --conf examples/s2m2_uncoupled.json --beta 123.45`
* `julia bin/sos.jl --conf examples/s2m2_uncoupled.json --beta 123.45 --basis-size 30`
* `julia bin/trotter.jl --conf examples/s2m2_uncoupled.json --beta 123.45 --basis-size 30 --num-links 200`
* `julia bin/sampling.jl --conf examples/s2m2_uncoupled.json --beta 123.45 --num-links 200 --num-samples 100000`


## Acknowledgements

Thanks to Neil Raymond for designing the parameter file format and helping to verify this implementation!


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
