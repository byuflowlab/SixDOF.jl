# Reference frames

Besides the six-degree-of-freedom integrator, SixDOF.jl holds the rigid-body
kinematic tree that the FLOW Lab's aerodynamic solvers share:
[LiftingLines.jl](https://github.com/byuflowlab/LiftingLines.jl),
[VortexLattice.jl](https://github.com/byuflowlab/VortexLattice.jl) and
[FLOWPanel.jl](https://github.com/byuflowlab/FLOWPanel.jl). Each of those used
to carry its own copy of the same tree; having one here is what lets a
lifting-line rotor, a vortex-lattice wing and a panelled fuselage sit in one
simulation and move together. The tree prescribes motion; it does not solve
for it. Coupling it to the integrator, so the vehicle's state comes from the
bodies' forces, is the next step.

A solver's user rarely calls the tree directly: LiftingLines' `simulate!`
takes the frames and moves the bodies itself. The tree is called directly when
writing a maneuver, or when writing a new body type.

## Building a tree

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

The examples below use the smallest body that can join the tree: one marked
point with a hub and an axis. (The protocol it implements is in the next
section; `examples/frames.jl` is the runnable version of all of this.)

```julia
using SixDOF, StaticArrays, LinearAlgebra

struct Marker{TF}
    p::SVector{3,TF}      # the point
    hub::SVector{3,TF}    # where it spins about
    axis::SVector{3,TF}   # and about which axis
end
Base.eltype(::Marker{TF}) where TF = TF
Base.eltype(::Type{Marker{TF}}) where TF = TF
SixDOF.move(m::Marker, origin, R, dx) =
    Marker(R * (m.p - origin) + origin + dx, R * (m.hub - origin) + origin + dx, R * m.axis)
SixDOF.spin(m::Marker) = (m.hub, m.axis)
```

A vehicle flying at 20 m/s along `-x` with two rotors of 0.5 m radius at
`y = ±1 m`, spinning about `x` at 100 rad/s in opposite senses. The root owns
nothing directly; each rotor gets a child frame at its hub:

```julia
bodies = Any[Marker(SVector(0.0,  1.5, 0.0), SVector(0.0,  1.0, 0.0), SVector(1.0, 0.0, 0.0)),
             Marker(SVector(0.0, -1.5, 0.0), SVector(0.0, -1.0, 0.0), SVector(1.0, 0.0, 0.0))]

frames = ReferenceFrame(bodies; v = [-20.0, 0.0, 0.0], dependent_index = Int[])
add_frame!(frames, "left",  1, [0.0,  1.0, 0.0], [1]; omega_axis = [1.0, 0.0, 0.0], omega =  100.0)
add_frame!(frames, "right", 1, [0.0, -1.0, 0.0], [2]; omega_axis = [1.0, 0.0, 0.0], omega = -100.0)

[f.name for f in frames]          # ["vehicle", "left", "right"]
frames[1].dependent_index          # [1, 2]: the root carries both rotors along
frames[2].dependent_index          # [1]
```

`rotor_frames` builds exactly this from the bodies' own `spin`:

```julia
frames = rotor_frames(bodies; omega = [100.0, -100.0], v = [-20.0, 0.0, 0.0])
[f.name for f in frames]          # ["vehicle", "rotor_1", "rotor_2"]
```

## Reading the motion

What a solver needs from the tree each step is how fast each body is moving,
to impose the kinematic velocity on the flow.

```@docs
frame_motion
point_velocity
```

For the left rotor, the chain gives the hub's position, the flight speed it
inherits from the vehicle, and its own shaft rate:

```julia
origin, v, omega = frame_motion(frames, 1)
origin                             # [0.0, 1.0, 0.0]
v                                  # [-20.0, 0.0, 0.0]
omega                              # [100.0, 0.0, 0.0]

u = point_velocity(frames, 1, bodies[1].p)
u                                  # [-20.0, 0.0, 50.0]: flight speed plus 100 rad/s x 0.5 m
norm(u)                            # 53.85
```

The fluid velocity the point sees on account of its own motion is `-u`; that
is what a solver subtracts at its control points.

## Moving the bodies

```@docs
propagate_kinematics!
```

Call it at the end of a time step, after the loads have been taken and the
wake shed at the current position. The integration is explicit Euler over
`dt` in each frame. Nine steps of a quarter revolution each 1/36 of a turn:

```julia
dt = 2pi / 100 / 36                # 10 degrees of rotor per step
for k in 1:9
    propagate_kinematics!(bodies, frames, dt)
end
bodies[1].p                        # [-0.314, 1.0, 0.5]: a quarter turn about the hub, 0.314 m of flight
frames[1].x                        # [-0.314, 0.0, 0.0]
frames[2].R                        # the left frame's basis, turned 90 degrees about x
```

The tree moves a body through a two-function protocol that the body's own
package implements by qualified name; the tree never looks inside a body.

```@docs
SixDOF.move
SixDOF.spin
```

A body must also define `Base.eltype` (its number type) so that
`ReferenceFrame(bodies)` can pick the tree's number type; `ReferenceFrame(TF, n)`
builds the root without bodies. LiftingLines' `Line` and `Rotor` implement the
protocol by rewriting their section arrays in place and returning a fresh
immutable wrapper, so a device buffer that holds the sections is never
reallocated.

## Maneuvers

Frames are immutable; a maneuver replaces a frame with a copy whose motion has
changed, between steps:

```@docs
ReferenceFrame(::ReferenceFrame)
```

A pitch-up of 10 degrees per second beginning at `t = 0.05 s`, applied before
each step:

```julia
function pitch_up!(frames, t)
    frames[1] = ReferenceFrame(frames[1]; omega_axis = [0.0, 1.0, 0.0], omega = t > 0.05 ? deg2rad(10) : 0.0)
    return nothing
end

t = 0.0
for k in 1:20
    pitch_up!(frames, t)
    propagate_kinematics!(bodies, frames, 0.01)
    t += 0.01
end
frames[1].R * [1.0, 0.0, 0.0]      # [0.9997, 0.0, -0.0244]: the nose is 1.4 degrees up
frame_motion(frames, 1)[3]         # [99.97, 0.17, -2.44]: the rotor's rate now has the pitch rate in it
```

A solver's `simulate!` takes a `maneuver!` callback and calls it with the
frames each step, so this is the whole of what a maneuver is.

## Rotations

```@docs
Rodrigues
inverse_Rodrigues
```
