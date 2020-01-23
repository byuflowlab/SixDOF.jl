module SixDOF

using LinearAlgebra: norm, cross

export Control, MassProp, Reference
export AbstractAeroModel, AbstractPropulsionModel, AbstractInertialModel, 
    AbstractAtmosphereModel, AbstractController
export StabilityDeriv, MotorPropBatteryDataFit, UniformGravitationalField, 
    ConstantAtmosphere, ConstantController
export CO, COUNTER, COCOUNTER
export sixdof!


# ------ General Structs -------

"""
    State(x, y, z, phi, theta, psi, u, v, w, p, q, r)

State of the aircraft: positions in inertial frame, euler angles,
velocities in body frame, angular velocities in body frame.
"""
struct State{TF}
    x::TF  # position (inertial frame)
    y::TF
    z::TF
    phi::TF  # orientation, euler angles
    theta::TF
    psi::TF
    u::TF  # velocity (body frame)
    v::TF
    w::TF
    p::TF  # angular velocity (body frame)
    q::TF
    r::TF
end

"""
    Control(de, dr, da, df, throttle)

Define the control settings: delta elevator, delta rudder, delta aileron, 
delta flaps, and throttle.
"""
struct Control{TF}
    de::TF  # elevator
    dr::TF  # rudder
    da::TF  # aileron
    df::TF  # rudder
    throttle::TF
end

"""
    MassProp(m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz)

Mass and moments of inertia in the body frame.
Ixx = int(y^2 + z^2, dm)
Ixz = int(xz, dm)

Most aircraft are symmetric about y and so there is a convenience method 
to specify only the nonzero components.
MassProp(m, Ixx, Iyy, Izz, Ixz)
"""
struct MassProp{TF}
    m::TF
    Ixx::TF 
    Iyy::TF
    Izz::TF
    Ixz::TF 
    Ixy::TF
    Iyz::TF
end

# most aircraft are symmetric in y
MassProp(m, Ixx, Iyy, Izz, Ixz) = MassProp(m, Ixx, Iyy, Izz, Ixz, zero(Ixx), zero(Ixx))  


"""
    Reference(S, b, c)

The reference area, span, and chord used in the aerodynamic computations.
"""
struct Reference{TF}
    S::TF  # area
    b::TF  # span
    c::TF  # chord
end

# ----------------------------------------------


# --------------- Interfaces ---------------

abstract type AbstractAtmosphereModel end

"""
    wind(model::AbstractAtmosphereModel, state)

Compute wind velocities.

**Returns**
- Wi: wind velocities in inertial frame
- Wb: gust velocities in body frame (just a convenience to allow some velocities in body frame)
"""
function wind(model::AbstractAtmosphereModel, state)
    @warn "wind function not implemented for AbstractAtmosphereModel"
    Wi = [0.0, 0.0, 0.0]
    Wb = [0.0, 0.0, 0.0]
    return Wi, Wb
end

"""
    properties(model::AbstractAtmosphereModel, state)

Compute atmospheric density and the speed of sound.
"""
function properties(model::AbstractAtmosphereModel, state)
    @warn "properties function not implemented for AbstractAtmosphereModel"
    rho = 1.225  # sea-level properties
    asound = 340.3
    return rho, asound
end

"""
    gravity(model::AbstractAtmosphereModel, state)

Compute the local acceleration of gravity.
"""
function gravity(model::AbstractAtmosphereModel, state)
    @warn "gravity function not implemented for AbstractAtmosphereModel"
    g = 9.81
    return g
end


# ----

abstract type AbstractAeroModel end

"""
    aeroforces(model::AbstractAeroModel, atm::AbstractAtmosphereModel, state::State, control::Control, mp::MassProp, ref::Reference)

Compute the aerodynamic forces and moments in the body frame.
return F, M
"""
function aeroforces(model::AbstractAeroModel, atm, state, control, mp, ref)
    @warn "aeroforces function not implemented for AbstractAeroModel"
    # forces and moments in body frame
    F = [0.0, 0.0, 0.0]
    M = [0.0, 0.0, 0.0]
    return F, M
end

# ----

abstract type AbstractPropulsionModel end

"""
    propulsionforces(model::AbstractPropulsionModel, atm::AbstractAtmosphereModel, state::State, control::Control, mp::MassProp, ref::Reference)

Compute the propulsive forces and moments in the body frame.
return F, M
"""
function propulsionforces(model::AbstractPropulsionModel, atm, state, control, mp, ref)
    @warn "propulsionforces function not implemented for AbstractPropulsionModel"
    # forces and moments in body frame
    F = [0.0, 0.0, 0.0]
    M = [0.0, 0.0, 0.0]
    return F, M
end

# ----

abstract type AbstractInertialModel end

"""
    gravityforces(model::AbstractInertialModel, atm::AbstractAtmosphereModel, state::State, control::Control, mp::MassProp, ref::Reference)

Compute the gravitational forces and moments in the body frame.
return F, M
"""
function gravityforces(model::AbstractInertialModel, atm, state, control, mp, ref)
    @warn "gravityforces function not implemented for AbstractInertialModel"
    # forces and moments in body frame
    F = [0.0, 0.0, 0.0]
    M = [0.0, 0.0, 0.0]
    return F, M
end


# ----

abstract type AbstractController end

"""
    setcontrol(controller::AbstractController, time, atm::AbstractAtmosphereModel, state::State, lastcontrol::Control, mp::MassProp, ref::Reference)

Compute control state for next time step given current state  
return control::Control
"""
function setcontrol(controller::AbstractController, time, atm, state, mp, ref)
    @warn "setcontrol function not implemented for AbstractController"
    control = Control(0.0, 0.0, 0.0, 0.0, 0.0)
    return control
end


# -----------------------------




# ------------- helper functions (private) --------------


"""
    inertialtobody(state)

Construct a rotation matrix from inertial frame to body frame

The assumed order of rotation is 
1) psi radians about the z axis, 
2) theta radians about the y axis, 
3) phi radians about the x axis. 

This is an orthogonal transformation so its inverse is its transpose.
"""
function inertialtobody(state)

    R = Array{eltype(state.phi)}(undef, 3, 3)

    cphi, ctht, cpsi = cos.([state.phi, state.theta, state.psi])
    sphi, stht, spsi = sin.([state.phi, state.theta, state.psi])

    R[1, 1] = ctht*cpsi
    R[1, 2] = ctht*spsi
    R[1, 3] = -stht

    R[2, 1] = sphi*stht*cpsi - cphi*spsi
    R[2, 2] = sphi*stht*spsi + cphi*cpsi
    R[2, 3] = sphi*ctht

    R[3, 1] = cphi*stht*cpsi + sphi*spsi
    R[3, 2] = cphi*stht*spsi - sphi*cpsi
    R[3, 3] = cphi*ctht

    return R

end


"""
    windtobody(alpha, beta)

Rotation matrix from wind frame to body frame.
- alpha: angle of attack
- beta: sideslip angle
"""
function windtobody(alpha, beta)

    ca, cb = cos.([alpha, beta])
    sa, sb = sin.([alpha, beta])

    Rwb = [ca*cb  -ca*sb  -sa;
           sb      cb     0.0;
           sa*cb  -sa*sb  ca]

    return Rwb
end

"""
    stabilitytobody(alpha, beta)

Rotation matrix from stability frame to body frame.
- alpha: angle of attack
"""
function stabilitytobody(alpha)

    ca = cos(alpha)
    sa = sin(alpha)

    Rsb = [ca  0.0  -sa;
           0.0 1.0  0.0;
           sa  0.0  ca]

    return Rsb
end


"""
    windaxes(atm::AbstractAtmosphereModel, state)

Compute relative velocity in wind axes (airspeed, aoa, sideslip)
"""
function windaxes(atm, state)

    # velocity vectors
    Vb = [state.u, state.v, state.w]
    Wi, Wb = wind(atm, state)

    Rib = inertialtobody(state)

    # relative wind
    Vrel = Vb - (Rib*Wi + Wb)

    # airspeed
    Va = norm(Vrel)

    # angle of attack
    alpha = atan(Vrel[3], Vrel[1])

    # sideslip
    beta = asin(Vrel[2] / Va)

    return Va, alpha, beta

end


# ----------------------------------------------------


# ------- Some Default Interface Implementations -----

"""
    StabilityDeriv(CL0, CLalpha, CLq, CLM, CLdf, CLde, alphas, 
        CD0, U0, exp_Re, e, Mcc, CDdf, CDde, CDda, CDdr, 
        CYbeta, CYp, CYr, CYda, CYdr, Clbeta, 
        Clp, Clr, Clda, Cldr, 
        Cm0, Cmalpha, Cmq, CmM, Cmdf, Cmde, 
        Cnbeta, Cnp, Cnr, Cnda, Cndr)

Stability derivatives of the aircraft.  Most are self explanatory if you are
familiar with stability derivatives (e.g., CLalpha is dCL/dalpha or the lift curve slope).
Some less familiar ones include
- M: Mach number
- alphas: the angle of attack for stall
- U0: the speed for the reference Reynolds number CD0 was computed at
- exp_Re: the coefficient in the denominator of the skin friction coefficient (0.5 laminar, 0.2 turbulent)
- e: Oswald efficiency factor
- Mcc: crest critical Mach number (when compressibility drag rise starts)

"""
struct StabilityDeriv{TF} <: AbstractAeroModel
    CL0::TF
    CLalpha::TF
    CLq::TF
    CLM::TF
    CLdf::TF
    CLde::TF
    alphas::TF  # TODO: should probably do in terms of CLmax

    CD0::TF
    U0::TF  # velocity corresponding to Reynolds number of CD0  (TODO: rethink this)
    exp_Re::TF  # exponent for Reynolds number scaling. typical values: exp_Re = 0.5 laminar, 0.2 turbulent
    e::TF  # Oswald efficiency factor
    Mcc::TF  # crest-critical Mach number when compressibility drag rise starts (quartic)
    CDdf::TF
    CDde::TF
    CDda::TF
    CDdr::TF

    CYbeta::TF
    CYp::TF
    CYr::TF
    CYda::TF
    CYdr::TF

    Clbeta::TF
    Clp::TF
    Clr::TF
    Clda::TF
    Cldr::TF

    Cm0::TF
    Cmalpha::TF
    Cmq::TF
    CmM::TF
    Cmdf::TF
    Cmde::TF

    Cnbeta::TF
    Cnp::TF
    Cnr::TF
    Cnda::TF
    Cndr::TF
end


"""
A simple (mostly) linear aerodynamics model
"""
function aeroforces(sd::StabilityDeriv, atm, state, control, ref, mp)

    # airspeed, angle of attack, sideslip
    Va, alpha, beta = windaxes(atm, state)

    # Mach number and dynamic pressure
    rho, asound = properties(atm, state)
    Mach = Va / asound
    qdyn = 0.5 * rho * Va^2

    # rename for convenience
    p = state.p
    q = state.q
    r = state.r
    de = control.de
    df = control.df
    dr = control.dr
    da = control.da


    # lift
    CL = sd.CL0 + sd.CLalpha*alpha + sd.CLq*q *ref.c/(2*Va) + sd.CLM*Mach
        + sd.CLdf*df + sd.CLde*de
    
    em = exp(-50*(alpha - sd.alphas))
    ep = exp(50*(alpha + sd.alphas))
    sigma = (1 + em + ep)/((1 + em)*(1 + ep))
    CL = (1- sigma)*CL + sigma * 2 * sign(alpha)*sin(alpha)^2*cos(alpha)

    # drag
    CDp = sd.CD0*(Va/sd.U0)^sd.exp_Re  
    CDi = CL^2/(pi*(ref.b^2/ref.S)*sd.e)
    CDc = (Mach < sd.Mcc) ? 0.0 : 20*(Mach - sd.Mcc)^4 

    CD = CDp + CDi + CDc + abs(sd.CDdf*df) + abs(sd.CDde*de) + abs(sd.CDda*da) + abs(sd.CDdr*dr)

    # side force
    CY = sd.CYbeta*beta + (sd.CYp*p + sd.CYr*r)*ref.b/(2*Va) + sd.CYda*da + sd.CYdr*dr
    
    # rolling moment
    Cl = sd.Clbeta*beta + (sd.Clp*p + sd.Clr*r)*ref.b/(2*Va) + sd.Clda*da + sd.Cldr*dr

    # pitching moment
    Cm = sd.Cm0 + sd.Cmalpha*alpha + sd.Cmq*q * ref.c/(2*Va) + sd.CmM*Mach + sd.Cmdf*df + sd.Cmde*de

    # yawing moment
    Cn = sd.Cnbeta*beta + (sd.Cnp*p + sd.Cnr*r)*ref.b/(2*Va) + sd.Cnda*da + sd.Cndr*dr

    # transfer forces from wind to body axes
    Rwb = windtobody(alpha, beta)

    F = Rwb*[-CD, CY, -CL] * qdyn * ref.S

    M = Rwb*[Cl*ref.b, Cm*ref.c, Cn*ref.b] * qdyn * ref.S

    return F, M
end

@enum PropType CO=1 COUNTER=-1 COCOUNTER=0

"""
    MotorPropBatteryDataFit(CT2, CT1, CT0, CQ2, CQ1, CQ0, D, num, type,
        R, Kv, i0, voltage)

**Inputs**
- CT2, CT1, CT0: quadratic fit to propeller thrust coefficient of form: CT = CT2*J2 + CT1*J + CT0
- CQ2, CQ1, CQ0: quadratic fit to propeller torque coefficient of form: CQ = CQ2*J2 + CQ1*J + CQ0
- D: propeller diameter
- num: number of propellers
- type: CO (torques add), COUNTER (torques add but with minus sign), COCOUNTER (no torque, they cancel out)
- R: motor resistance
- Kv: motor Kv
- i0: motor no-load current
- voltage: battery voltage
"""
struct MotorPropBatteryDataFit{TF, TI, PropType} <: AbstractPropulsionModel
    # CT = CT2*J2 + CT1*J + CT0
    # CQ = CQ2*J2 + CQ1*J + CQ0
    CT2::TF  # prop data fit
    CT1::TF
    CT0::TF
    CQ2::TF
    CQ1::TF
    CQ0::TF
    D::TF  # prop diameter
    num::TI
    type::PropType

    R::TF  # motor resistance
    Kv::TF  # motor Kv
    i0::TF  # motor no-load current

    voltage::TF  # battery voltage
end

function propulsionforces(prop::MotorPropBatteryDataFit, atm, state, control, ref, mp)

    # airspeed, angle of attack, sideslip
    Va, _, _ = windaxes(atm, state)

    # density
    rho, _ = properties(atm, state)

    D = prop.D

    # determine torque for motor/prop match (quadratic equation)
    a = rho*D^5/(2*pi)^2 * prop.CQ0
    b = rho*D^4/(2*pi)*Va * prop.CQ1 + 1.0/(prop.R*prop.Kv)
    c = rho*D^3*Va^2 * prop.CQ2 - control.throttle*prop.voltage/(prop.R*prop.Kv) + prop.i0/prop.Kv
    Omega = (-b + sqrt(b^2 - 4*a*c))/(2*a)

    # advance ratio
    n = Omega/(2*pi)
    J = Va/(n*D)

    # thrust and torque
    CT = prop.CT0 + prop.CT1*J + prop.CT2*J^2
    CQ = prop.CQ0 + prop.CQ1*J + prop.CQ2*J^2

    T = prop.num * CT * rho * n^2 * D^4
    Q = prop.num * CQ * rho * n^2 * D^5 * Int(prop.type)

    return [T, 0, 0], [Q, 0, 0] 
end

"""
    UniformGravitationalField()

Assumes center of mass and center of gravity are coincident.
"""
struct UniformGravitationalField <: AbstractInertialModel end

function gravityforces(model::UniformGravitationalField, atm, state, control, ref, mp)

    W = mp.m * gravity(atm, state)
    ct, cp = cos.([state.theta, state.phi])
    st, sp = sin.([state.theta, state.phi])

    Fg = W*[-st, ct*sp, ct*cp]
    Mg = [zero(W), zero(W), zero(W)]  # no gravitational moment

    return Fg, Mg
end


"""
    ConstantAtmosphere(Wi, Wb, rho, asound, g)

Constant atmospheric properties.
"""
struct ConstantAtmosphere{TF, TV<:AbstractVector{TF}} <: AbstractAtmosphereModel
    Wi::TV
    Wb::TV
    rho::TF
    asound::TF
    g::TF
end


function wind(atm::ConstantAtmosphere, state)
    return atm.Wi, atm.Wb
end

function properties(atm::ConstantAtmosphere, state)
    return atm.rho, atm.asound
end

function gravity(atm::ConstantAtmosphere, state)
    return atm.g
end


"""
    ConstantController(de, dr, da, df, throttle)

Just a dummy controller that outputs constant control outputs the whole time.
"""
struct ConstantController{TF} <: AbstractController
    de::TF
    dr::TF
    da::TF
    df::TF
    throttle::TF
end

function setcontrol(controller::ConstantController, time, atm, state, mp, ref)
    return Control(controller.de, controller.dr, controller.da, controller.df, controller.throttle)
end
# --------------------------------------------------------


# ------------- main functions (public) --------------

"""
    sixdof!(ds, s, params, time)

dynamic and kinematic ODEs.  Follows format used in DifferentialEquations package.
- s = x, y, z, phi, theta, psi, u, v, w, p, q, r (same order as State)
- params = control, massproperties, reference, aeromodel, propmodel, inertialmodel, atmmodel
"""
function sixdof!(ds, s, params, time)

    x, y, z, phi, theta, psi, u, v, w, p, q, r = s
    mp, ref, aeromodel, propmodel, inertialmodel, atmmodel, controller = params

    # ---- controller -------
    state = State(s...)
    control = setcontrol(controller, time, atmmodel, state, mp, ref)
    # -----------------------
    
    # --------- forces and moments ---------
    # aerodynamics
    Fa, Ma = aeroforces(aeromodel, atmmodel, state, control, ref, mp)

    # propulsion
    Fp, Mp = propulsionforces(propmodel, atmmodel, state, control, ref, mp)

    # weight
    Fg, Mg = gravityforces(inertialmodel, atmmodel, state, control, ref, mp)

    # total forces and moments
    F = Fa + Fp + Fg
    M = Ma + Mp + Mg

    # --------------------------------------


    # ----- derivative of state --------
    Vb = [u, v, w]
    omegab = [p, q, r]

    # linear kinematics
    Rib = inertialtobody(state)
    rdot = Rib' * Vb

    # angular kinematics
    phidot = p + (q*sin(phi) + r*cos(phi))*tan(theta)
    thetadot = q*cos(phi) - r*sin(phi)
    psidot = (q*sin(phi) + r*cos(phi))/cos(theta)

    # linear dynamics
    vdot = F/mp.m - cross(omegab, Vb)

    # angular dynamics
    I = [mp.Ixx -mp.Ixy -mp.Ixz;
         -mp.Iyz mp.Iyy -mp.Iyz;
         -mp.Ixz -mp.Iyz mp.Izz]
    omegadot = I \ (M - cross(omegab, I*omegab))

    # -------------------------

    # TODO: if we need more efficiency we can avoid allocating and then assigning.
    ds[1:3] = rdot
    ds[4] = phidot
    ds[5] = thetadot
    ds[6] = psidot
    ds[7:9] = vdot
    ds[10:12] = omegadot
end


end # module
