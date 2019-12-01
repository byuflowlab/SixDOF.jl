var documenterSearchIndex = {"docs":
[{"location":"theory/#Theory-1","page":"Theory","title":"Theory","text":"","category":"section"},{"location":"theory/#Introduction-1","page":"Theory","title":"Introduction","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"This nonlinear sixDOF model uses aircraft conventions as that is our target application.  Although the model is fairly general, some extensions may be needed for other applications (derivatives of moment of inertia tensor for morphing aircraft, quaternions and gravitational moments for spacecraft).","category":"page"},{"location":"theory/#State-and-Coordinate-Systems-1","page":"Theory","title":"State and Coordinate Systems","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Two main coordinate systems are used in this analysis as depicted below.  The first is the inertial frame (x_i y_i z_i), fixed to the ground [1], centered at some arbitrary reference point (sometimes this coordinate system is referred to as north-east-down).  The second is the body frame (x_b y_b z_b), fixed to the moving and rotating vehicle, centered at the vehicle's center of mass.  A third coordinate system will be introduced later in connection with aerodynamics.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"(Image: coordinate systems)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The aircraft has 12 state variables: the linear and angular positions and the linear and angular velocities.  ","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"s = x y z phi theta psi u v w p q r^T","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The positions are defined in the inertial frame, while the velocities are defined in the body frame.  The rotational velocities p q r correspond to roll, pitch, and yaw rates respectively.  ","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Although there are twelve state variables, there are only six indepedent degrees of freedom (hence the name SixDOF).  The other six constraints are purely kinematic relationships enforcing consistency between positions and velocities.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"We define the angular orientation of the vehicle using Euler angles.  Euler angles are not unique. We use a typical representation with the following order from the inertial axis to the body axis:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"rotate psi radians about the inertial z-axis\nrotate theta radians about the new y-axis\nrotate phi radians about the new x-axis.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The resulting rotations will then align with the body axis.  The two coordinate systems will generally still differ by a translation.  For the calculations here, like transfering velocity vectors between the frames, the translation offset is irrelevant.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The resulting rotation matrix from inertial to body is (see Beard and McLain for derivation [2]:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"R_i rightarrow b = \nbeginbmatrix\ncosthetacospsi  costhetasinpsi  -sintheta \nsinphisinthetacospsi - cosphisinpsi  sinphisinthetasinpsi + cosphicospsi  sinphicostheta \ncosphisinthetacospsi + sinphisinpsi  cosphisinthetasinpsi - sinphicospsi  cosphicostheta\nendbmatrix","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"This matrix is orthogonal so its inverse is its transpose.  In other words the rotation matrix from body to inertial is the transpose of above:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"R_b rightarrow i = R_i rightarrow b^T","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The downside of Euler angles is that there is a singularity.  For this particular order that singularity occurs at a theta (associated with pitch) of pm 90^circ.  At that state, the yaw angle is not unique, a phenomenon known as gimbal lock.  The avoid this issue, quaternions are commonly used.  For spacecraft such pitch angles would be normal and thus this would be an important consideration.  For now, we are only simulating aircraft so the simplicity of Euler angles is justified.","category":"page"},{"location":"theory/#Kinematics-and-Dynamics-1","page":"Theory","title":"Kinematics and Dynamics","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"We now need to compute derivatives of all state variables in order to formulate an ODE.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The first two relationships are simply kinematic relationships between positions and velocities.  The derivative of the vehicle's position (vecr_i = x y z^T) is just the vehicle's velocity (vecV_b = u v w^T).  However, position is defined in the inertial frame, whereas the velocity is given in the body frame.  Thus, a rotation matrix is needed.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"fracd vecr_idt = R_b rightarrow i vecV_b","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Next, we need the derivatives of the Euler angles.  The procedure is the same as above, except that the Euler angles are defined across four different coordinate systems consisting of the body frame, inertial frame, and two intermediate frames between those two (see above definition of Euler angles).  The details are just coordinate transformations, but are a bit tedious and so are not repeated here (see [2] for derivation).  The result is:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"beginbmatrix\ndotphi \ndottheta \ndotpsi \nendbmatrix\n=\nbeginbmatrix\n1  sinphitantheta  cosphitantheta \n0  cosphi  -sinphi \n0  sinphicostheta  cosphicostheta \nendbmatrix\nbeginbmatrix\np \nq \nr \nendbmatrix","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Now we need to resolve the dynamics. We start with the linear velocity: vecV_b = u v w^T.  We can apply Newton's second law to the motion of the vehicle as follows:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"vecF = m fracd vecVdt_i","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"where the derivative is in an inertial frame.  Newton's laws are only applicable in an inertial frame of reference. If a reference frame b is rotating relative to an inertial frame i then the time derivative of any vector vecv can be described as:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"leftfracd vecvdtright_i + leftfracd vecvdtright_b + vecomega_bi times vecv","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"where vecomega_bi is the rotation of reference frame b relative to i.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"We transform Newton's law to the body frame using the relationship from above.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"vecF = m left( fracdvecVdt_b + vecomega_bi times vecV right)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"We can evaluate these forces/velocities in any frame, but the most convenient will be the body frame:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"vecF_b = m left( fracdvecV_bdt_b + vecomega_b times vecV_b right)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Since the end goal is an ODE and so we solve for dV_bdt_b:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"fracd vecV_bdt_b = fracvecF_bm - vecomega_b times  vecV_b","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"where vecomega_b = p q r^T.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"We follow a similar procedure for the rotational motion.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"beginaligned\nvecM = fracdvecHdt_i\n = fracdvecHdt_b + vecomega_bi times vecH \n = fracdvecH_bdt_b + vecomega_b times vecH_b \nendaligned","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The angular momentum is given by:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"beginaligned\nvecH = I vecomega\nvecH_b = I_b vecomega_b\nendaligned","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"where I_b is the moment of inertia tensor defined in the body axis.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"I = \nbeginbmatrix\nI_xx  -I_xy  -I_xz \n-I_xy  I_yy  -I_yz \n-I_zx  -I_zy  I_zz \nendbmatrix","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"(note the minus signs on the off-diagonal components).  The diagonal components of the tensor are given by:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"I_xx = int (y^2 + z^2) dm","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"and the off-diagonal components are given by:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"I_xy = int xy dm","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Most aircraft are symmetric about the y-axis and so typically I_xy = I_yz = 0.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"We now plug the definition of the angular moment back into Newton's law.  Unless, the aircraft is able to morph, the inertia tensor is constant in the body frame and so the equation becomes:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"vecM_b = I_b fracdvecomega_bdt_b + vecomega_b times I_b vecomega_b","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"We solve for the derivative of the angular velocity:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"fracdvecomega_bdt_b = I_b^-1 left(vecM_b - vecomega_b times I_b vecomega_b right)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Note, that although a matrix inverse was written, in numerical implementation a linear solve would be used rather than an explicit inversion.","category":"page"},{"location":"theory/#Forces/Moments-1","page":"Theory","title":"Forces/Moments","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The primary forces and moments on the aircraft are from aerodynamics, propulsion, and gravity.  The former two are outside the scope of the dynamics module, although simple default implementations will be discussed.  The latter is straightforward.  Gravity always acts about the vehicle's center of mass [3], and so cannot create any moments on the vehicle.  Based on our above derivation we need forces/moments to be specified in the body frame.","category":"page"},{"location":"theory/#Gravity-1","page":"Theory","title":"Gravity","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"In the inertial frame gravity is in the positive z direction (recall that positive z is down).  Thus we just need to apply a rotation to the body frame:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"F_g = R_i rightarrow b\nbeginbmatrix\n0 \n0 \nmg \nendbmatrix","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"For our choice of Euler angles we can write this out explicitly as:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"beginbmatrix\n-mg sintheta \nmg sinphicostheta \nmg cosphicostheta \nendbmatrix","category":"page"},{"location":"theory/#Aerodynamics-(General)-1","page":"Theory","title":"Aerodynamics (General)","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Aerodynamic behavior is relative to the local freestream V_infty.  The local freestream, is the negative of what is referred to as true airspeed, and the true airspeed is ground speed minus wind speed.  ","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"vecV_infty = -vecV_a = -(vecV_g - vecV_w)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Ground speed, in the body frame is vecV_g = vecV_b = u v w^T.  For convenience, we allow the user to define a portion of the wind in the inertial frame and a portion in the body frame if desired.  The former would typically be steady-state winds, and the latter would typically represent gusts.  While an aerodynamicist uses the freestream, we will use dynamics convention and use the vehicle's airspeed instead (just a change in sign):","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"vecV_a_b =  vecV_b - left( R_i rightarrow b vecV_wi + vecV_wbright)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Aerodynamics requires an additional coordinate system: the wind frame.  Actually, we need one more intermediate frame, the stability frame, which occurs between the body frame and the wind frame.  Again, we adopt dynamics convention (aerodynamics convention uses the opposite sign for the x and z axes).  These two coordinate systems are depicted below.  The top figure shows the stability frame.  It is rotated about the body y-axis axis by the angle of attack alpha.  The bottom two figures depict the wind axis.  It is rotated about the stability z-axis by the sideslip angle beta.  The coordinate transformation from wind to body axes (the direction we will usually go) is thus:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"beginbmatrix\ncosalpha cosbeta  -cosalpha sinbeta  -sinalpha \nsinbeta  cosbeta  0 \nsinalpha cosbeta  -sinalpha sinbeta  cosalpha \nendbmatrix","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"(and the transformation from body to wind axes is the transpose of this).","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"(Image: Stability and wind axes)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The airspeed vector in the wind axes is:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"vecV_a_w = (V_a 0 0)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"From the definitions, either using the figure, or the coordinate transformation we can see that the airspeed vector in the body frame is then:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"vecV_a_b = V_a \nbeginbmatrix\ncosalphacosbeta\nsinbeta \nsinalphacosbeta\nendbmatrix","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"We see that the airspeed is the magnitude of this vector","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"V_a = vecV_a_b_2","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"If we call the components of the airspeed in the body frame as follows: vecV_a_b = u_a v_a w_w^T then the angle of attack and sideslip angles are given by:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"alpha = tan^-1left(fracw_au_aright)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"beta = sin^-1left(fracv_aV_aright)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"These calculations may also be visualized from the figure below.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"(Image: angle of attack and sideslip)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Because we have used dynamic conventions, we need to be careful to translate the aerodynamic forces properly.  Lift, side force, and drag, would be defined with the following signs using the dynamics wind frame:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"F_aero_w = -D Y -L","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Notice that drag and lift are defined in opposite directions from the dynamics positive coordinate directions.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Aerodynamic moments are typically also defined in the stability or wind axes, and so would be defined as follows using the dynamics wind frame.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"M_aero_w = mathcalL M N","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"where the moments are rolling, pitching, and yaw respectively. Note that the rolling moment and yawing moment do not require any sign changes as the aerodynamic directions follow the same as used in dynamics.","category":"page"},{"location":"theory/#Aerodynamics-(a-specific-implementation)-1","page":"Theory","title":"Aerodynamics (a specific implementation)","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Everything up to this point defines the 6DOF simulator.  The aerodynamics, propulsion, and other forces/moments if applicable, are specific to the model choice.  As a simple default a conventional ``linear'' aerodynamics model is included.  This model we use is actually not strictly linear, particuarly in the case of drag, which is fundamentally nonlinear.  Arguably, that addition is negliglbe as aircraft drag typically plays a negligible role in its dynamic behavior.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"We define the aerodynamic forces and moments using stability derivatives and control derivatives. The control deflections include delta f delta e delta r delta a for flaps, elevator, rudder, and aileron respectively.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"For lift:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"C_L = C_L_0 + fracd C_Ldalphaalpha + fracd C_Ldq fracq c2 V_a + fracd C_LdMM + fracd C_Ld delta fdelta f + fracd C_Ld delta edelta e","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"where Mach is the Mach number.  A maximum and minimum C_L is also enforced (a very crude approximation of stall) via C_L_max and C_L_min.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Drag is a combination of parasitic, induced, and compressibility drag.  The parasitic drag coefficient is proportional to a skin friction coefficient C_D_f, which is Reynolds number dependent.  For laminar flow (Blasisus solution):","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"C_D_f propto frac1sqrtRe","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"For turbulent flow no such analytic solution exists, various empirical relationships exists such as the Schlichting formula:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"C_D_f propto frac1Re^02","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"We assume that the kinematic viscosity does not change appreciably throughout the trajectory and so the Reynolds number scaling only changes with the flight speed.  We can then define the parasitic drag as:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"C_D_p = C_D_0 left(fracV_aV_refright)^-k","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"where C_D_0 is the zero-lift drag coefficient, and the second term is the Reynolds number correction based on the reference speed V_ref (corresponding to the Reynolds number at which C_D_0 was computed at, and an exponent for the skin friction coefficnet: k = 05 for laminar and k = 02 for turbulent.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The induced drag is given by its standard definition:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"C_D_i = fracC_L^2pi AR e","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"where e is Oswald's efficiency factor, and the aspect ratio uses the reference dimensions of the aircraft (span and area).","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Compressibility drag uses an simple empirical quartic drag rise after the crest-critical Mach number:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"C_D_c = 20(M - M_cc)^4","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The parasitic drag is the sum of these three components plus that from control deflections.  An absolute value must be used for the control parts as deflections in either direction would increase drag.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"C_D = C_D_p + C_D_i + C_D_c + leftfracd C_Dd delta fdelta fright + leftfracd C_Dd delta edelta eright + leftfracd C_Dd delta adelta aright + leftfracd C_Dd delta rdelta rright","category":"page"},{"location":"theory/#Propulsion-1","page":"Theory","title":"Propulsion","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"[1]: The earth is not, strictly speaking, an inertial frame since it is rotating and so objects on the surface are accelerating.  However, for our applications including the inertial effect is negligible.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"[2]: Small Unmanned Aircraft: Theory and Practice, Randal W. Beard and Timothy W. McLain, Princeton University Press, 2012.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"[3]: Actually it acts about the center of gravity, but in practice these are the same.  The only time we'd see a noticeable difference is if the vehicle was extremely large, or if this was a satellite then the very small gravitational moments might be important in its dynamics.","category":"page"},{"location":"#SixDOF.jl-1","page":"Home","title":"SixDOF.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [SixDOF]","category":"page"},{"location":"#SixDOF.aeroforces-Tuple{StabilityDeriv,AtmosphereModel,Any,Any,Any}","page":"Home","title":"SixDOF.aeroforces","text":"A simple (mostly) linear aerodynamics model\n\n\n\n\n\n","category":"method"},{"location":"#SixDOF.inertialtobody-Tuple{Any}","page":"Home","title":"SixDOF.inertialtobody","text":"Construct a rotation matrix from inertial frame to body frame\n\nThe assumed order of rotation is \n\npsi radians about the z axis, \ntheta radians about the y axis, \nphi radians about the x axis. \n\nThis is an orthogonal transformation so its inverse is its transpose.\n\n\n\n\n\n","category":"method"},{"location":"#SixDOF.wind-Tuple{AtmosphereModel,Any}","page":"Home","title":"SixDOF.wind","text":"Returns\n\nu, v, w: wind velocities in inertial frame\nug, vg, wg: gust velocities in body frame (just a convenience to allow some velocities in body frame)\n\n\n\n\n\n","category":"method"},{"location":"#SixDOF.windaxes-Tuple{AtmosphereModel,Any}","page":"Home","title":"SixDOF.windaxes","text":"Compute relative velocity in wind axes (airspeed, aoa, sideslip)\n\n\n\n\n\n","category":"method"}]
}
