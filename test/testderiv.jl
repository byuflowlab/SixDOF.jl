using SixDOF
import ForwardDiff
import ReverseDiff
# using FLOWMath: centraldiff


function wrapper(x)
    num = 2
    type = COCOUNTER

    mp = MassProp(x[1:5]...)
    ref = Reference(x[6:8]...)
    sd = StabilityDeriv(x[9:45]...)
    propulsion = MotorPropBatteryDataFit(x[46], x[47], x[48], x[49], x[50], x[51], x[52], num, type, x[53:56]...)
    atm = ConstantAtmosphere(x[57:59], x[60:62], x[63], x[64], x[65])
    controller = ConstantController(x[66:70]...)
    inertial = UniformGravitationalField()
    Vinf = 10.0
    alpha = 3.0*pi/180
    s0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Vinf*cos(alpha), 0.0, Vinf*sin(alpha), 0.0, 0.0, 0.0]
    time = 1.0
    params = mp, ref, sd, propulsion, inertial, atm, controller
    ds = zeros(eltype(x[1]), 12)
    sixdof!(ds, s0, params, time)
    return ds
end

x = zeros(71)

x[1] = 1.56
x[2] = 0.1147
x[3] = 0.0576
x[4] = 0.1712
x[5] = 0.0015
x[6] = 0.2589
x[7] = 1.4224
x[8] = 0.3302
x[9] = 0.09167 # Zero-alpha lift
x[10] = 3.5016  # lift curve slope
x[11] = 2.8932 # Pitch rate derivative
x[12] = 0.0 # Mach derivative
x[13] = 0.0  # flaps derivative
x[14] = 0.2724  # elevator derivative
x[15] = 20*pi/180
x[16] = 0.01631  # zero-lift drag coerff
x[17] = 10.0  # velocity corresponding to Reynolds number of CD0
x[18] = -0.2  # exponent in Reynolds number calc
x[19] = 0.8  # Oswald efficiency
x[20] = 0.7  # crest critcal Mach number
x[21] = 0.0  # flaps
x[22] = 0.3045  # elevator
x[23] = 0.0  # aileron
x[24] = 0.0  # rudder
x[25] = -0.07359 # Sideslip derivative
x[26] = 0.0  # Roll rate derivative
x[27] = 0.0 # Yaw rate derivative
x[28] = 0.0 # Roll control (aileron) derivative
x[29] = 0.0 # Yaw control (rudder) derivative
x[30] = -0.02854  # Sideslip derivative
x[31] = -0.3209  # Roll rate derivative
x[32] = 0.03066  # Yaw rate derivative
x[33] = 0.1682  # Roll (aileron) control derivative
x[34] = 0.0  #Yaw (rudder) control derivative
x[35] = -0.02338 # Zero-alpha pitch
x[36] = -0.5675 # Alpha derivative
x[37] = -1.3990 # Pitch rate derivative
x[38] = 0.0
x[39] = 0.0
x[40] = -0.3254 # Pitch control derivative
x[41] = -0.00040  # Slideslip derivative
x[42] = -0.01297  # Roll rate derivative
x[43] = -0.00434  # Yaw rate derivative
x[44] = -0.00328  # Roll (aileron) control derivative
x[45] = 0.0  # Yaw (rudder) control derivative
x[46] = 0.11221
x[47] = -0.13803
x[48] = -0.047394
x[49] = 0.0062
x[50] = 0.00314
x[51] = -0.015729
x[52] = 10*0.0254
x[53] = 0.5
x[54] = 2500.0 * pi/30
x[55] = 0.3
x[56] = 8.0
x[57:59] .= [0.0, 0.0, 0.0]
x[60:62] .= [0.0, 0.0, 0.0]
x[63] = 1.2682
x[64] = 300.0
x[65] = 9.81
x[66:70] .= [-3.0*pi/180, 0.0, 0.0, 0.0, 0.8]

f = wrapper(x)

J = ForwardDiff.jacobian(wrapper, x)
J2 = ReverseDiff.jacobian(wrapper, x)
# J3 = centraldiff(wrapper, x)

# diff = J - J3
# idx = (diff .!= 0.0)
# rerr = diff[idx]./J3[idx]

# @test maximum(abs.(rerr)) < 1e-5

# diff = J2 - J3
# idx = (diff .!= 0.0)
# rerr = diff[idx]./J3[idx]

# @test maximum(abs.(rerr)) < 1e-5

diff = J - J2
idx = (diff .!= 0.0)
rerr = diff[idx]./J[idx]
# println(maximum(abs.(rerr)))
@test maximum(abs.(rerr)) < 1e-7

# @code_warntype wrapper(x)