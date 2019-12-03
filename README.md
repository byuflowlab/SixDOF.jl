# SixDOF.jl

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://andrewning.github.io/SixDOF.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/SixDOF.jl)

*A nonlinear six degree of freedom dynamics model*
Author: Andrew Ning

![coordinate systems](docs/src/cs.png)

The dynamics implementation is quite standard and could be used for any rigid body moving through 3D space, although our usage (and conventions) are for aircraft.  The main limitation of the current model is that it is based on Euler angles, rather than quaternions, but that again is because of our applications.  All forces/moments and the atmosphere model are abstract types to allow for user-defined models.  The default implementations of these abstract types include:
- aerodynamics: a standard aircdraft stability/control-derivative based model (the drag model is more detailed, but that generally has a negligible impact on the dynamics).
- propulsion: an electric motor/propeller/battery model using a first-order motor model and nondimensional propeller thrust and torque curves.
- gravity: center of mass and center of gravity are coincident.
- atmosphere: constant properties across the trajectory.

Author: Andrew Ning

### Install

```julia
pkg> add https://github.com/byuflowlab/SixDOF.jl
```

### Run Unit Tests

```julia
pkg> activate .
pkg> test
```

### Getting Started

Users should read the Guide section of the [documentation](https://flow.byu.edu/SixDOF.jl).  Developers may be interested in the Theory section of the documentation.