# Guide

## Basic Usage

```@setup zagi
using PyPlot
```

First, we import the module.

```@example zagi
using SixDOF
```

There are several structs we need to define.  The inputs used in this example correspond to the Zagi flying wing in Appendix E of Small Unmanned Aircraft: Theory and Practice by Beard and McLain.  We specify the mass properties:

```@docs
MassProp
```

For this example:
```@example zagi
m = 1.56
Ixx = 0.1147
Iyy = 0.0576
Izz = 0.1712
Ixz = 0.0015
mp = MassProp(m, Ixx, Iyy, Izz, Ixz)
nothing # hide
```

Next, we specify reference areas and lengths:
```@docs
Reference
```

```@example zagi
Sref = 0.2589
bref = 1.4224
cref = 0.3302
ref = Reference(Sref, bref, cref)
nothing # hide
```

A controller builds off the abstract type `AbstractController`, which must supply the following method:
```@docs
SixDOF.setcontrol
```

There is a default implementation in the module, which is just a dummy controller that outputs contant control output
```@docs
ConstantController
```

In this example, we don't use any control deflections.  Just throttle, at 80%.
```@example zagi
controller = ConstantController(0.0, 0.0, 0.0, 0.0, 0.8)
nothing # hide
```

Next, we define the atmospheric model.  `AbstractAtmosphereModel` is an abstract type that must define the following three methods.
```@docs
SixDOF.wind(::AbstractAtmosphereModel, ::Any)
SixDOF.properties(::AbstractAtmosphereModel, state)
SixDOF.gravity
```
There is a default implementation in the module, which is the simplest possible model: one with constant properties:
```@docs
ConstantAtmosphere
```

We use that in this example:
```@example zagi
Wi = [0.0, 0.0, 0.0]
Wb = [0.0, 0.0, 0.0]
rho = 1.2682
asound = 300.0
g = 9.81
atm = ConstantAtmosphere(Wi, Wb, rho, asound, g)
nothing # hide
```

Finally, we now need to define the forces and moments.  We provide three abstract types for three types of forces: aerodynamics, propulsion, and gravity.  In principle you could use these methods to define forces/moments for any application, but for aircraft this is a convenient separation.

All three forces take in all the same inputs, which include everything discussed so far and the state.  The state is an internally used struct that contains the state of the aircraft (or other object).
```@docs
SixDOF.State
```

The `AbstractAeroModel` abstract type must define the following function:
```@docs
SixDOF.aeroforces(model::AbstractAeroModel, atm, state, control, mp, ref)
```
The default implementation of `AbstractAeroModel` is one based on stability derivatives.
```@docs
StabilityDeriv
```

We use the following values for this example.
```@example zagi

CL0 = 0.09167 # Zero-alpha lift
CLalpha = 3.5016  # lift curve slope
CLq = 2.8932 # Pitch rate derivative
CLM = 0.0 # Mach derivative
CLdf = 0.0  # flaps derivative
CLde = 0.2724  # elevator derivative
CLmax = 1.4  # max CL (stall)
CLmin = -0.9  # min CL (neg stall)
alphas = 20*pi/180

CD0 = 0.01631  # zero-lift drag coerff
U0 = 10.0  # velocity corresponding to Reynolds number of CD0
exp_Re = -0.2  # exponent in Reynolds number calc
e = 0.8  # Oswald efficiency
Mcc = 0.7  # crest critcal Mach number
CDdf = 0.0  # flaps
CDde = 0.3045  # elevator
CDda = 0.0  # aileron
CDdr = 0.0  # rudder

CYbeta = -0.07359 # Sideslip derivative
CYp = 0.0  # Roll rate derivative
CYr = 0.0 # Yaw rate derivative
CYda = 0.0 # Roll control (aileron) derivative
CYdr = 0.0 # Yaw control (rudder) derivative

Clbeta = -0.02854  # Sideslip derivative
Clp = -0.3209  # Roll rate derivative
Clr = 0.03066  # Yaw rate derivative
Clda = 0.1682  # Roll (aileron) control derivative
Cldr = 0.0  #Yaw (rudder) control derivative

Cm0 = -0.02338 # Zero-alpha pitch
Cmalpha = -0.5675 # Alpha derivative
Cmq = -1.3990 # Pitch rate derivative
CmM = 0.0
Cmdf = 0.0
Cmde = -0.3254 # Pitch control derivative

Cnbeta = -0.00040  # Slideslip derivative
Cnp = -0.01297  # Roll rate derivative
Cnr = -0.00434  # Yaw rate derivative
Cnda = -0.00328  # Roll (aileron) control derivative
Cndr = 0.0  # Yaw (rudder) control derivative

sd = StabilityDeriv(CL0, CLalpha, CLq, CLM, CLdf, CLde, alphas,
    CD0, U0, exp_Re, e, Mcc, CDdf, CDde, CDda, CDdr,
    CYbeta, CYp, CYr, CYda, CYdr,
    Clbeta, Clp, Clr, Clda, Cldr,
    Cm0, Cmalpha, Cmq, CmM, Cmdf, Cmde,
    Cnbeta, Cnp, Cnr, Cnda, Cndr)
nothing # hide
```


The `AbstractPropulsionModel` abstract type must define the following function:
```@docs
SixDOF.propulsionforces(model::AbstractPropulsionModel, atm, state, control, mp, ref)
```
The default implementation of `AbstractPropulsionModel` is based on a first-order motor model coupled with a parameterized curve fit of propeller data.  The torque for when the motor and propeller are matched is solved for and then used to compute thrust.
```@docs
MotorPropBatteryDataFit
```


This example uses data roughly corresponding to an APC thin electric 10x5
```@example zagi

CT0 = 0.11221
CT1 = -0.13803
CT2 = -0.047394
CQ0 = 0.0062
CQ1 = 0.00314
CQ2 = -0.015729
D = 10*0.0254
num = 2
type = COCOUNTER
R = 0.5
Kv = 2500.0 * pi/30
i0 = 0.3
voltage = 8.0
propulsion = MotorPropBatteryDataFit(CT2, CT1, CT0, CQ2, CQ1, CQ0, D, num, type, R, Kv, i0, voltage)
nothing # hide
```

Finally, the `AbstractInertialModel` must implment the following function
```@docs
SixDOF.gravityforces(model::AbstractInertialModel, atm, state, control, mp, ref)
```

The default implementation (``UniformGravitationalField``) assumes that the center of mass and center of gravity are coincident and so there is no gravitational moment.  The default will likely be used most of the time as that condition is true for almost all applications, except perhaps some spacecraft in high orbits where small gravitational torques may matter.

```@example zagi
inertial = UniformGravitationalField()
nothing # hide
```

The main function is 
```@docs
sixdof!
```

We rarely use it directly, but rather use it in connection with an ODE solver.  In this case the ``DifferentialEquations`` package.  We start with an initial velocity at an angle of attack and simulate for four seconds.

```@example zagi
import DifferentialEquations

Vinf = U0
alpha = 3.0*pi/180
s0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Vinf*cos(alpha), 0.0, Vinf*sin(alpha), 0.0, 0.0, 0.0]
tspan = (0.0, 4.0)
p = mp, ref, sd, propulsion, inertial, atm, controller


prob = DifferentialEquations.ODEProblem(sixdof!, s0, tspan, p)
sol = DifferentialEquations.solve(prob)
nothing # hide
```

We can plot the results.  For example the linear positions and velocities.  The y-components are not plotted in this case, because they are all zero as there are no control deflections or wind that would cause lateral motion in this example.
```@example zagi
using PyPlot

figure()
plot(sol.t, sol[1, :])
xlabel("time (s)")
ylabel("x inertial position (m)")
savefig("x.svg"); nothing # hide
figure()
plot(sol.t, sol[3, :])
xlabel("time (s)")
ylabel("z inertial position (m)")
savefig("z.svg"); nothing # hide
figure()
plot(sol.t, sol[7, :])
xlabel("time (s)")
ylabel("u body velocity (m/s)")
savefig("u.svg"); nothing # hide
figure()
figure()
plot(sol.t, sol[9, :])
xlabel("time (s)")
ylabel("w body velocity (m/s)")
savefig("w.svg"); nothing # hide
```
![](x.svg)
![](z.svg)
![](u.svg)
![](w.svg)

## Additional Helper Functions

A few other helper functions exist.  These first three create rotation matrices from one coordinate system to the body frame.  All of these rotation matrices are orthogonal so the inverse of the matrix is just its transpose.

```@docs
SixDOF.inertialtobody
SixDOF.windtobody
SixDOF.stabilitytobody
```

Another function extracts airspeed, angle of attack, and sideslip from the body motion and wind.

```@docs
SixDOF.windaxes
```

