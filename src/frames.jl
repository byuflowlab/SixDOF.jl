#=
Rigid-body kinematic tree, shared by the aerodynamic solvers (LiftingLines,
VortexLattice, FLOWPanel): a flat `Vector{ReferenceFrame}` whose root is the
vehicle, children are rotors, wings, anything that moves relative to its
parent. The tree knows nothing about the bodies it moves: each body type
implements the two-function protocol below.

    SixDOF.move(body, origin, R, dx) -> body      rotate about `origin` by `R`, then translate by `dx`
    SixDOF.spin(body) -> nothing | (origin, axis)  where a body spins in place, for `rotor_frames`

and `eltype(body)` gives its number type. The two are deliberately not
exported: a solver extends them by qualified name and never calls them itself. A solver takes the motion of a body
from [`frame_motion`](@ref) or the velocity of one of its points from
[`point_velocity`](@ref); a body moving at `U` sees fluid velocity `-U`.
=#

"""
    SixDOF.move(body, origin, R, dx)

Rigid motion of a body: rotate it about the global point `origin` by the
rotation matrix `R`, then translate it by `dx`; return the moved body (which may
be the same object mutated, or a new immutable wrapper). Every body type the
frame tree moves adds a method, `function SixDOF.move(b::MyBody, origin, R, dx)`;
[`propagate_kinematics!`](@ref) is the only caller. Not exported.
"""
function move end

"""
    SixDOF.spin(body)

`nothing` for a body that does not spin in place, or `(origin, axis)` in global
coordinates for one that does (a rotor: its hub and shaft axis). Used by
[`rotor_frames`](@ref); the default is `nothing`, so only spinning bodies add
a method. Not exported.
"""
spin(::Any) = nothing


"""
    Rodrigues(axis, angle)

Rotation matrix for a right-handed rotation of `angle` [rad] about the unit
vector `axis`.
"""
function Rodrigues(axis, angle::TF) where TF
    s, c = sincos(angle)
    t = one(TF) - c
    x, y, z = axis
    return SMatrix{3,3,TF,9}(
        t*x*x + c,   t*x*y + s*z, t*x*z - s*y,
        t*x*y - s*z, t*y*y + c,   t*y*z + s*x,
        t*x*z + s*y, t*y*z - s*x, t*z*z + c,
    )
end

"""
    inverse_Rodrigues(R)

Recover the axis-angle vector (axis scaled by the angle) of a rotation matrix.
Returns the zero vector for the identity.
"""
function inverse_Rodrigues(R::SMatrix{3,3,TF,9}) where TF
    arg = (R[1,1] + R[2,2] + R[3,3] - one(TF)) / 2
    theta = acos(clamp(arg, -one(TF), one(TF)))
    s = sin(theta)
    iszero(s) && return zero(SVector{3,TF})
    return SVector{3,TF}(R[3,2] - R[2,3], R[1,3] - R[3,1], R[2,1] - R[1,2]) *
           (theta / (2 * s))
end

#--- the frame ---#

"""
    ReferenceFrame{TF}

One node of a rigid-body kinematic tree. Frames are stored flat, in a
`Vector{ReferenceFrame}` whose first entry is the root; parent and child links
are indices into that vector.

**Fields** -- `x`, `v`, `omega_axis` and `R` are all expressed in the *parent*
frame; only `Rp2g` reaches global coordinates.
- `x`: origin
- `v`: translational velocity
- `omega_axis`, `omega`: rotation axis (unit) and rate [rad/s]
- `R`: this frame's basis
- `Rp2g`: the parent's basis in global coordinates, cached by
  [`propagate_kinematics!`](@ref)
- `name`: used by the `String` forms of [`add_frame!`](@ref) and
  [`frame_index`](@ref)
- `parent_index`: `-1` for the root
- `child_index`: child frames
- `dependent_index`: geometries attached to this frame; a parent's list is the
  union of its descendants', maintained by `add_frame!`
"""
struct ReferenceFrame{TF}
    x::SVector{3,TF}
    v::SVector{3,TF}
    omega_axis::SVector{3,TF}
    omega::TF
    R::SMatrix{3,3,TF,9}
    Rp2g::SMatrix{3,3,TF,9}
    name::String
    parent_index::Int
    child_index::Vector{Int}
    dependent_index::Vector{Int}
end

Base.eltype(::ReferenceFrame{TF}) where TF = TF
Base.eltype(::Type{ReferenceFrame{TF}}) where TF = TF

function ReferenceFrame{TF}(f::ReferenceFrame) where TF
    return ReferenceFrame{TF}(
        SVector{3,TF}(f.x), SVector{3,TF}(f.v), SVector{3,TF}(f.omega_axis),
        TF(f.omega), SMatrix{3,3,TF,9}(f.R), SMatrix{3,3,TF,9}(f.Rp2g),
        f.name, f.parent_index, copy(f.child_index), copy(f.dependent_index),
    )
end
Base.convert(::Type{ReferenceFrame{TF}}, f::ReferenceFrame) where TF =
    ReferenceFrame{TF}(f)

"""
    ReferenceFrame(f::ReferenceFrame; x, v, omega_axis, omega, R)

A copy of `f` with the given fields replaced (all in the parent frame). This is
how a maneuver changes a frame's motion between steps, since frames are
immutable: `frames[1] = ReferenceFrame(frames[1]; omega = 0.2)`.
"""
function ReferenceFrame(f::ReferenceFrame{TF}; x = f.x, v = f.v, omega_axis = f.omega_axis,
                        omega = f.omega, R = f.R) where TF
    return ReferenceFrame{TF}(SVector{3,TF}(x), SVector{3,TF}(v), _unit(SVector{3,TF}(omega_axis)),
                              TF(omega), SMatrix{3,3,TF,9}(R), f.Rp2g, f.name, f.parent_index,
                              f.child_index, f.dependent_index)
end

"""
    ReferenceFrame(geometries; kwargs...)

Build the root frame of a kinematic tree and return it as a one-element
`Vector`, the container every other function here takes. Attach children with
[`add_frame!`](@ref).

`geometries` is the `Vector` of bodies the tree moves (only their `eltype` and
count are read); by default the root owns all of them. `ReferenceFrame(TF, n;
kwargs...)` builds the same root for `n` bodies of number type `TF`.

# Keyword Arguments
- `origin`, `v`: origin and velocity in global coordinates (default zero)
- `omega_axis`, `omega`: rotation axis and rate (default `+z`, zero)
- `R`, `Rp2g`: bases (default identity)
- `name`: default `"vehicle"`
- `dependent_index`: geometries owned by the root
"""
function ReferenceFrame(geometries::AbstractVector;
        origin = nothing,
        v = nothing,
        omega_axis = nothing,
        omega = nothing,
        R = nothing,
        Rp2g = nothing,
        name::String = "vehicle",
        child_index::Vector{Int} = Int[],
        dependent_index::Vector{Int} = collect(eachindex(geometries)),
    )
    TF = isempty(geometries) ? Float64 : promote_type(eltype.(geometries)...)
    return ReferenceFrame(TF, length(geometries); origin, v, omega_axis, omega, R, Rp2g,
                          name, child_index, dependent_index)
end

function ReferenceFrame(::Type{TF0}, n::Integer;
        origin = nothing,
        v = nothing,
        omega_axis = nothing,
        omega = nothing,
        R = nothing,
        Rp2g = nothing,
        name::String = "vehicle",
        child_index::Vector{Int} = Int[],
        dependent_index::Vector{Int} = collect(1:n),
    ) where TF0
    TF = promote_type(TF0, _eltype_or(origin, TF0), _eltype_or(v, TF0),
                      _eltype_or(omega_axis, TF0), _eltype_or(omega, TF0),
                      _eltype_or(R, TF0), _eltype_or(Rp2g, TF0))
    return [ReferenceFrame{TF}(
        isnothing(origin) ? zero(SVector{3,TF}) : SVector{3,TF}(origin),
        isnothing(v) ? zero(SVector{3,TF}) : SVector{3,TF}(v),
        isnothing(omega_axis) ? SVector{3,TF}(0, 0, 1) : _unit(SVector{3,TF}(omega_axis)),
        isnothing(omega) ? zero(TF) : TF(omega),
        isnothing(R) ? _identity(TF) : SMatrix{3,3,TF,9}(R),
        isnothing(Rp2g) ? _identity(TF) : SMatrix{3,3,TF,9}(Rp2g),
        name, -1, child_index, dependent_index,
    )]
end

_eltype_or(::Nothing, ::Type{TF}) where TF = TF
_eltype_or(x, ::Type{TF}) where TF = eltype(x)
_identity(::Type{TF}) where TF = SMatrix{3,3,TF,9}(1, 0, 0, 0, 1, 0, 0, 0, 1)
_unit(v::SVector{3,TF}) where TF = (n = norm(v); iszero(n) ? v : v / n)

"""
    rotor_frames(geometries; omega, kwargs...)

The frame tree for rotors that each spin in place: a root frame (with the
`kwargs` of [`ReferenceFrame`](@ref), e.g. a vehicle velocity) and one child
per body that [`spin`](@ref)s, at its own hub rotating about its own axis at
`omega` (one rate, or one per geometry). Bodies that do not spin stay on the root. Building
`ReferenceFrame(geometries; omega_axis, omega)` instead spins every geometry
about the root origin, which is what one rotor at the origin wants and what
several rotors side by side do not.
"""
function rotor_frames(geometries::AbstractVector; omega, kwargs...)
    frames = ReferenceFrame(geometries; dependent_index = Int[], kwargs...)
    om = omega isa AbstractVector ? omega : fill(omega, length(geometries))
    length(om) == length(geometries) || throw(DimensionMismatch("one omega per geometry is required"))
    root_dependents = Int[]
    for (i, g) in enumerate(geometries)
        sp = spin(g)
        if sp !== nothing
            add_frame!(frames, "rotor_$i", 1, sp[1], [i]; omega_axis = sp[2], omega = om[i])
        else
            push!(root_dependents, i)
        end
    end
    isempty(root_dependents) || _update_dependent_indices!(frames, 1, root_dependents)
    return frames
end

"""
    frame_motion(frames, geometry_index)

The rigid motion of the frame that owns geometry `geometry_index`, in global
coordinates: `(origin, v, omega_vec)` such that a point `p` of the geometry
moves at `v + omega_vec x (p - origin)`. Accumulated down the tree, so a rotor
on a translating wing frame carries the wing's velocity and its own shaft rate.
Throws if no frame owns the geometry.
"""
function frame_motion(frames::Vector{<:ReferenceFrame{TF}}, geometry_index::Int) where TF
    return _frame_motion(frames, 1, geometry_index, zero(SVector{3,TF}), _identity(TF),
                         zero(SVector{3,TF}), zero(SVector{3,TF}))
end

function _frame_motion(frames, i_frame, gi, dx_p2g::SVector{3,TF}, R_p2g, v_parent, omega_parent) where TF
    frame = frames[i_frame]
    origin = R_p2g * frame.x + dx_p2g
    v = v_parent + cross(omega_parent, origin - dx_p2g) + R_p2g * frame.v
    omega_vec = omega_parent + (R_p2g * frame.omega_axis) * frame.omega
    gi in frame.dependent_index && !any(c -> gi in _owned_below(frames, c), frame.child_index) &&
        return (origin, v, omega_vec)
    R = R_p2g * frame.R
    for c in frame.child_index
        gi in _owned_below(frames, c) || continue
        return _frame_motion(frames, c, gi, origin, R, v, omega_vec)
    end
    i_frame == 1 && throw(ArgumentError("no frame owns geometry $gi"))
    return (origin, v, omega_vec)
end

# Geometry indices owned by frame `i` or any frame below it.
function _owned_below(frames, i)
    owned = copy(frames[i].dependent_index)
    for c in frames[i].child_index
        append!(owned, _owned_below(frames, c))
    end
    return owned
end

"""
    frame_index(frames, name)

Index of the frame called `name`. Throws if there is no such frame.
"""
function frame_index(frames::AbstractVector{<:ReferenceFrame}, name::AbstractString)
    for (i, frame) in enumerate(frames)
        frame.name == name && return i
    end
    throw(ArgumentError("no reference frame named \"$name\""))
end

"""
    add_frame!(frames, name, parent, origin, geometry_indices; kwargs...)

Append a child frame and return its index. `parent` is either the parent's index
or its name, and `origin` is the child's origin *in the parent frame*. The
geometries listed are registered as dependents of the new frame and of every
ancestor, so that moving an ancestor carries them along.

Keyword arguments are those of [`ReferenceFrame`](@ref), all expressed in the
parent frame.
"""
function add_frame!(frames::Vector{ReferenceFrame{TF}}, name::String,
        parent_index::Int, origin, geometry_indices::Vector{Int};
        v = zero(SVector{3,TF}),
        omega_axis = SVector{3,TF}(0, 0, 1),
        omega = zero(TF),
        R = _identity(TF),
        Rp2g = _identity(TF),
    ) where TF
    1 <= parent_index <= length(frames) ||
        throw(ArgumentError("parent_index $parent_index is out of range"))
    push!(frames, ReferenceFrame{TF}(
        SVector{3,TF}(origin), SVector{3,TF}(v), _unit(SVector{3,TF}(omega_axis)),
        TF(omega), SMatrix{3,3,TF,9}(R), SMatrix{3,3,TF,9}(Rp2g),
        name, parent_index, Int[], copy(geometry_indices),
    ))
    push!(frames[parent_index].child_index, length(frames))
    _update_dependent_indices!(frames, parent_index, geometry_indices)
    return length(frames)
end

add_frame!(frames::Vector{<:ReferenceFrame}, name::String, parent_name::AbstractString,
           origin, geometry_indices::Vector{Int}; kwargs...) =
    add_frame!(frames, name, frame_index(frames, parent_name), origin,
               geometry_indices; kwargs...)

# Register `geometry_indices` as dependents of `frames[parent_index]` and of
# every ancestor, so moving an ancestor carries them along.
function _update_dependent_indices!(frames::Vector{<:ReferenceFrame},
        parent_index::Int, geometry_indices::Vector{Int})
    for i in geometry_indices
        i in frames[parent_index].dependent_index ||
            push!(frames[parent_index].dependent_index, i)
    end
    grandparent = frames[parent_index].parent_index
    grandparent == -1 ||
        _update_dependent_indices!(frames, grandparent, geometry_indices)
    return nothing
end

"""
    point_velocity(frames, geometry_index, p)

Velocity in global coordinates of the point `p` of geometry `geometry_index`,
`v + omega_vec x (p - origin)` from [`frame_motion`](@ref). The fluid velocity
the point sees on account of its own motion is the negative of this.
"""
function point_velocity(frames::Vector{<:ReferenceFrame}, geometry_index::Int, p)
    origin, v, omega_vec = frame_motion(frames, geometry_index)
    return v + cross(omega_vec, p - origin)
end

"""
    SixDOF.owns(frames, geometry_index)

Whether any frame of the tree owns geometry `geometry_index`. Not exported.
"""
owns(frames::Vector{<:ReferenceFrame}, geometry_index::Int) =
    !isempty(frames) && geometry_index in frames[1].dependent_index

#--- moving the geometry ---#

"""
    propagate_kinematics!(geometries, frames, dt)

Advance the kinematic tree by one step of `dt`: move every geometry to where its
frame chain puts it at `t + dt`, and update each frame's stored origin, basis and
cached parent-to-global rotation.

Call this at the *end* of a step, after the loads at the current position have
been taken and after the wake has been shed, so that the next step begins with
the geometry and its frames consistent.
"""
function propagate_kinematics!(geometries::AbstractVector,
        frames::Vector{<:ReferenceFrame{TF}}, dt) where TF
    isempty(frames) && return nothing
    _propagate_kinematics!(geometries, 1, frames, zero(SVector{3,TF}),
                           _identity(TF), TF(dt))
    return nothing
end

function _propagate_kinematics!(geometries, i_frame::Int, frames,
        dx_parent_to_global::SVector{3,TF}, R_parent_to_global::SMatrix{3,3,TF},
        dt::TF) where TF
    frame = frames[i_frame]

    origin_global = R_parent_to_global * frame.x + dx_parent_to_global
    dx = frame.v * dt
    dx_global = R_parent_to_global * dx
    dtheta = frame.omega * dt
    R_step = Rodrigues(frame.omega_axis, dtheta)
    R_step_global = Rodrigues(R_parent_to_global * frame.omega_axis, dtheta)

    # Every dependent, including those owned by descendants: a child's geometry
    # is rigidly attached to this frame, so it takes this frame's step here and
    # the child's own additional step when the recursion reaches it.
    for i in frame.dependent_index
        geometries[i] = move(geometries[i], origin_global, R_step_global, dx_global)
    end

    frames[i_frame] = ReferenceFrame{TF}(
        frame.x + dx, frame.v, frame.omega_axis, frame.omega, R_step * frame.R,
        R_parent_to_global, frame.name, frame.parent_index,
        frame.child_index, frame.dependent_index,
    )

    dx_child = origin_global + dx_global
    R_child = R_parent_to_global * R_step * frame.R
    for i in frame.child_index
        _propagate_kinematics!(geometries, i, frames, dx_child, R_child, dt)
    end
    return nothing
end
