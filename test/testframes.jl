# The kinematic tree on a minimal body: a set of points with the move/spin protocol.
using StaticArrays, LinearAlgebra

struct Dots{TF}
    p::Vector{SVector{3,TF}}
    hub::SVector{3,TF}
    axis::SVector{3,TF}
    spins::Bool
end
Base.eltype(::Dots{TF}) where TF = TF
Base.eltype(::Type{Dots{TF}}) where TF = TF
function SixDOF.move(d::Dots{TF}, origin, R, dx) where TF
    for i in eachindex(d.p); d.p[i] = R * (d.p[i] - origin) + origin + dx; end
    Dots{TF}(d.p, SVector{3,TF}(R * (d.hub - origin) + origin + dx), SVector{3,TF}(R * d.axis), d.spins)
end
SixDOF.spin(d::Dots) = d.spins ? (d.hub, d.axis) : nothing

ex, ey, ez = SVector(1.0,0,0), SVector(0,1.0,0), SVector(0,0,1.0)

@testset "Rodrigues" begin
    R = Rodrigues(ez, pi/2)
    @test R * ex ≈ ey
    @test inverse_Rodrigues(R) ≈ ez * pi/2
    @test inverse_Rodrigues(Rodrigues(ex, 0.0)) == zero(SVector{3,Float64})
end

@testset "tree bookkeeping" begin
    bodies = Any[Dots([ex], zero(ex), ez, false), Dots([ey], ey, ex, true)]
    frames = ReferenceFrame(bodies)
    @test length(frames) == 1 && frames[1].dependent_index == [1, 2] && frames[1].parent_index == -1
    @test eltype(frames[1]) == Float64
    frames = ReferenceFrame(Float32, 3)
    @test eltype(frames[1]) == Float32 && frames[1].dependent_index == [1, 2, 3]
    frames = ReferenceFrame(bodies; dependent_index = Int[])
    i = add_frame!(frames, "arm", 1, ey, [2]; omega_axis = ex, omega = 2.0)
    @test i == 2 && frame_index(frames, "arm") == 2
    @test frames[1].dependent_index == [2] && frames[1].child_index == [2]
    j = add_frame!(frames, "tip", "arm", ex, [1])
    @test frames[2].dependent_index == [2, 1] && frames[1].dependent_index == [2, 1]
    @test_throws ArgumentError frame_index(frames, "none")
    @test owns(frames, 1) && !owns(ReferenceFrame(bodies; dependent_index = Int[]), 1)
end

@testset "rotor_frames" begin
    bodies = Any[Dots([ex], zero(ex), ez, false), Dots([ey], ey, ex, true), Dots([ez], 2ez, ey, true)]
    frames = rotor_frames(bodies; omega = [0.0, 3.0, 5.0], v = ex)
    @test length(frames) == 3
    @test frames[1].dependent_index == [2, 3, 1]      # rotors registered by add_frame!, the line last
    @test frames[2].x == ey && frames[2].omega_axis == ex && frames[2].omega == 3.0
    @test frames[3].x == 2ez && frames[3].omega == 5.0
    @test_throws DimensionMismatch rotor_frames(bodies; omega = [1.0, 2.0])
end

@testset "propagate and motion" begin
    # a rotor at hub (1,0,0) spinning about z on a vehicle translating at v = (0,0,1)
    rotor = Dots([SVector(1.0, 0.5, 0.0)], ex, ez, true)
    bodies = Any[rotor]
    frames = rotor_frames(bodies; omega = pi, v = ez)
    origin, v, w = frame_motion(frames, 1)
    @test origin ≈ ex && v ≈ ez && w ≈ pi * ez
    p = bodies[1].p[1]
    @test point_velocity(frames, 1, p) ≈ ez + cross(pi * ez, p - ex)
    propagate_kinematics!(bodies, frames, 0.5)          # quarter turn, up by 0.5
    @test bodies[1].p[1] ≈ ex + 0.5 * ez + Rodrigues(ez, pi/2) * SVector(0.0, 0.5, 0.0)
    @test bodies[1].hub ≈ SVector(1.0, 0.0, 0.5)
    @test frames[1].x ≈ 0.5 * ez
    @test frames[2].Rp2g ≈ frames[1].R              # cached parent basis
    # after the move the motion query follows the moved hub
    origin, v, w = frame_motion(frames, 1)
    @test origin ≈ SVector(1.0, 0.0, 0.5)
    # a nested chain: root spinning about z carries the child's own rate
    b2 = Any[Dots([ex], ex, ez, true)]
    fr = ReferenceFrame(b2; omega_axis = ez, omega = 2.0, dependent_index = Int[])
    add_frame!(fr, "r", 1, ex, [1]; omega_axis = ez, omega = 1.0)
    _, v, w = frame_motion(fr, 1)
    @test w ≈ 3ez && v ≈ cross(2ez, ex)
    @test_throws ArgumentError frame_motion(fr, 2)
end
