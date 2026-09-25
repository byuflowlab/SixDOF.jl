#=
The reference-frame tree on its own: a vehicle carrying two rotors, moved for
one rotor revolution, with the motion of each body read back from the tree.

    julia --project=. examples/frames.jl

A body joins the tree by giving it `Base.eltype` and `SixDOF.move`, and
`SixDOF.spin` if it spins in place; `Blob` below is the smallest such body, a
set of points with a hub and an axis. An aerodynamic solver's wing or rotor
does the same with its own type.
=#
using SixDOF, StaticArrays, LinearAlgebra, Printf

struct Blob{TF}
    points::Vector{SVector{3,TF}}
    hub::SVector{3,TF}
    axis::SVector{3,TF}
    spins::Bool
end
Base.eltype(::Blob{TF}) where TF = TF
Base.eltype(::Type{Blob{TF}}) where TF = TF
function SixDOF.move(b::Blob{TF}, origin, R, dx) where TF
    for i in eachindex(b.points)
        b.points[i] = R * (b.points[i] - origin) + origin + dx
    end
    return Blob{TF}(b.points, SVector{3,TF}(R * (b.hub - origin) + origin + dx), SVector{3,TF}(R * b.axis), b.spins)
end
SixDOF.spin(b::Blob) = b.spins ? (b.hub, b.axis) : nothing

# two rotors of radius 0.5 m, tips marked by one point each, at y = +-1 m,
# spinning about x; a fuselage point that does not spin
ex, ey, ez = SVector(1.0, 0, 0), SVector(0, 1.0, 0), SVector(0, 0, 1.0)
bodies = Any[Blob([SVector(0.0,  1.0, 0.5)], ey, ex, true),
             Blob([SVector(0.0, -1.0, 0.5)], -ey, ex, true),
             Blob([SVector(0.5,  0.0, 0.0)], zero(ex), ez, false)]

# vehicle flying at 20 m/s along -x, rotors at 100 rad/s (a counter-rotating pair)
frames = rotor_frames(bodies; omega = [100.0, -100.0, 0.0], v = SVector(-20.0, 0.0, 0.0))
println("frames: ", join((f.name for f in frames), ", "))

for (i, b) in enumerate(bodies)
    origin, v, w = frame_motion(frames, i)
    u = point_velocity(frames, i, b.points[1])
    @printf("body %d: origin (%5.1f %5.1f %5.1f) m, v (%5.1f %5.1f %5.1f) m/s, omega %6.1f rad/s, point speed %6.1f m/s\n",
            i, origin..., v..., norm(w), norm(u))
end

# one revolution of the rotors in 360 steps
dt = 2pi / 100 / 360
for _ in 1:360
    propagate_kinematics!(bodies, frames, dt)
end
@printf("after one revolution: vehicle at x = %.3f m (expected %.3f), rotor-1 tip back at (%.3f %.3f %.3f)\n",
        frames[1].x[1], -20 * 2pi / 100, bodies[1].points[1]...)
