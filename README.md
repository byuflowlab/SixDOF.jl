# SixDOF.jl

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://andrewning.github.io/SixDOF.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/SixDOF.jl/dev)
![](https://github.com/byuflowlab/SixDOF.jl/workflows/Run%20tests/badge.svg)

*A nonlinear six degree of freedom dynamics model, particularly for aircraft*

Author: Andrew Ning

<img src="docs/src/cs.png" width="500">

The dynamics implementation is quite standard and could be used for any rigid body moving through 3D space, although our usage (and conventions) are for aircraft.  The main limitation of the current model is that it is based on Euler angles, rather than quaternions, but that again is because of our applications.  All forces/moments and the atmosphere model are abstract types to allow for user-defined models.  The default implementations of these abstract types include:
- aerodynamics: a standard aircraft stability/control-derivative based model (the drag model is more detailed, but that generally has a negligible impact on the dynamics).
- propulsion: an electric motor/propeller/battery model using a first-order motor model and nondimensional propeller thrust and torque curves.
- gravity: center of mass and center of gravity are coincident (constant gravitational field across the body).
- atmosphere: constant properties throughout trajectory.

### Install

```julia
pkg> add SixDOF
```

### Run Unit Tests

```julia
pkg> test SixDOF
```

### Getting Started

Users should read the Guide section of the [documentation](https://flow.byu.edu/SixDOF.jl).  Developers may be interested in the Theory section of the documentation.
