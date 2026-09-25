# Reference frames

Besides the six-degree-of-freedom integrator, SixDOF.jl holds the rigid-body
kinematic tree that the FLOW Lab's aerodynamic solvers (LiftingLines.jl,
VortexLattice.jl, FLOWPanel.jl) share, so that a wing, a rotor and a fuselage
from different solvers can move together in one simulation. The tree is
independent of the integrator: it prescribes motion; it does not yet solve for
it.

## The tree

A simulation's frames are one flat `Vector{ReferenceFrame}`. The first entry is
the root (the vehicle); every other frame names its parent by index, and a
frame owns the *bodies* attached to it by their indices into the simulation's
vector of bodies. A frame's origin, velocity, rotation axis and basis are all
expressed in its parent's frame.

```@docs
ReferenceFrame
add_frame!
frame_index
rotor_frames
```

A vehicle flying at 30 m/s along `-x` with one rotor at `(1, 0, 0)` spinning at
100 rad/s about `x`:

```julia
frames = ReferenceFrame(bodies; v = [-30.0, 0.0, 0.0], dependent_index = Int[])   # the root owns nothing directly
add_frame!(frames, "rotor", 1, [1.0, 0.0, 0.0], [1]; omega_axis = [1.0, 0.0, 0.0], omega = 100.0)
```

`rotor_frames(bodies; omega)` builds the same thing for every body that
[`SixDOF.spin`](@ref)s.

## Moving the bodies

The tree moves a body through a two-function protocol that the body's own
package implements by qualified name; the tree never looks inside a body.

```@docs
SixDOF.move
SixDOF.spin
```

A body must also define `Base.eltype` (its number type) for
`ReferenceFrame(bodies)` to pick the tree's number type; `ReferenceFrame(TF, n)`
builds the root without bodies.

```@docs
propagate_kinematics!
```

Call it at the end of a time step, after loads have been taken and the wake
shed at the current position. The integration is explicit Euler over `dt` in
each frame; a solver that needs the position at an intermediate time (a
multi-stage integrator) calls it with a partial step and then back.

## The motion of a body

What a solver needs from the tree each step is how fast each of its points is
moving, to impose the kinematic velocity on the flow.

```@docs
frame_motion
point_velocity
```

The fluid velocity a point sees on account of its own motion is the negative
of `point_velocity`.

## Rotations

```@docs
Rodrigues
inverse_Rodrigues
```
