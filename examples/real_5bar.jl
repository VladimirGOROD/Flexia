using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff

using StaticArrays

const g = 9.81

bd01 = Body2D(10, 500,length = 0.50)
bd02 = Body2D(10, 500,length = 0.50)
bd03 = Body2D(10, 500,length = 0.50)

bd1 = Body2D(10, 500, length = 1.56)
bd2 = Body2D(10, 500, length = 1.56)
bd3 = Body2D(10, 500, length = 0.50)
bd4 = Body2D(10, 500, length = 1.56)
bd5 = Body2D(10, 500, length = 1.56)

bd1.forces[2] = (x) -> -bd1.mass * g
bd2.forces[2] = (x) -> -bd2.mass * g
bd3.forces[2] = (x) -> -bd3.mass * g
bd4.forces[2] = (x) -> -bd4.mass * g
bd5.forces[2] = (x) -> -bd5.mass * g
bd01.forces[2] = (x) -> -bd01.mass * g
bd02.forces[2] = (x) -> -bd02.mass * g
bd03.forces[2] = (x) -> -bd03.mass * g

jnt01 = FixedJoint(bd01)

jnt11 = HingeJoint(bd01, bd1)
jnt2 = HingeJoint(bd1, bd2)
jnt3 = HingeJoint(bd2, bd3)
jnt31 = HingeJoint(bd3, bd03)
jnt4 = HingeJoint(bd3, bd4)
jnt5 = HingeJoint(bd4, bd5)
jnt21 = HingeJoint(bd5, bd02)

jnt02 = FixedJoint(bd02)
jnt03 = FixedJoint(bd03)

tcp11 = TorsionalSpring(bd01, bd1, 1100.,0.0, 100.)
tcp2 = TorsionalSpring(bd1, bd2, 1100.,0.0, 100.)
tcp3 = TorsionalSpring(bd2, bd3, 1100.,0.0, 100.)
tcp31 = TorsionalSpring(bd03, bd3, 1100.,0.0, 100.)
tcp4 = TorsionalSpring(bd3, bd4, 1100.,0.0, 100.)
tcp5 = TorsionalSpring(bd4, bd5, 1100.,0.0, 100.)
tcp21 = TorsionalSpring(bd02, bd5, 1100.,0.0, 100.)

set_position_on_first_body!(jnt11, SA[0.5, 0])
set_position_on_second_body!(jnt11, SA[-1.56, 0.])

set_position_on_first_body!(jnt2, SA[1.56, 0])
set_position_on_second_body!(jnt2, SA[-1.56, 0])

set_position_on_first_body!(jnt3, SA[1.56, 0])
set_position_on_second_body!(jnt3, SA[-0.5, 0])

set_position_on_first_body!(jnt4, SA[0.5, 0])
set_position_on_second_body!(jnt4, SA[-1.56, 0])

set_position_on_first_body!(jnt5, SA[1.56, 0])
set_position_on_second_body!(jnt5, SA[-1.56, 0])

set_position_on_first_body!(jnt21, SA[1.56, 0])
set_position_on_second_body!(jnt21, SA[0.5, 0])

set_position_on_first_body!(jnt31, SA[0., 0])
set_position_on_second_body!(jnt31, SA[0.5, 0])

jnt01.pos = SA[-0.5, 0.0]
jnt02.pos = SA[1.3531+1+1+0.5+1, 0.0]
jnt03.pos = SA[1.3531+1, 1.56+1.56+0.78+1.]

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)
add!(sys, bd3)
add!(sys, bd4)
add!(sys, bd5)
add!(sys, bd01)
add!(sys, bd02)
add!(sys, bd03)

add!(sys, jnt01)
add!(sys, jnt02)
add!(sys, jnt03)

add!(sys, jnt11)
add!(sys, jnt2)
add!(sys, jnt21)
add!(sys, jnt3)
add!(sys, jnt31)
add!(sys, jnt4)
add!(sys, jnt5)

add!(sys, tcp11)
add!(sys, tcp2)
add!(sys, tcp21)
add!(sys, tcp3)
add!(sys, tcp31)
add!(sys, tcp4)
add!(sys, tcp5)

if (!assemble!(sys))
    println("Assembling failed!")
end

func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

bd1_x_ind, bd1_y_ind, bd1_t_ind = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)
bd3_x_ind, bd3_y_ind, bd3_t_ind = get_body_position_dofs(sys, bd3)
bd4_x_ind, bd4_y_ind, bd4_t_ind = get_body_position_dofs(sys, bd4)
bd5_x_ind, bd5_y_ind, bd5_t_ind = get_body_position_dofs(sys, bd5)

bd01_x_ind, bd01_y_ind, _ = get_body_position_dofs(sys, bd02)
bd02_x_ind, bd02_y_ind, _ = get_body_position_dofs(sys, bd02)
bd03_x_ind, bd03_y_ind, _ = get_body_position_dofs(sys, bd03)

initial = zeros(number_of_dofs(sys))

initial[bd01_x_ind] = 0
initial[bd01_y_ind] = 0

initial[bd1_x_ind] = 0
initial[bd1_y_ind] = 1.56
initial[bd1_t_ind] = pi/2

initial[bd2_x_ind] = 1.
initial[bd2_y_ind] = 1.56+1.56+0.78
initial[bd2_t_ind] = pi/4

initial[bd3_x_ind] = 1.3531+1
initial[bd3_y_ind] = 1.56+1.56+0.78+1

initial[bd4_x_ind] = 1.3531+1+1+0.5
initial[bd4_y_ind] = 2.6631+0.78+0.5
initial[bd4_t_ind] = - pi/4

initial[bd5_x_ind] = 1.3531+1+1+0.5+1
initial[bd5_y_ind] = 1.56
initial[bd5_t_ind] = -pi/2

initial[bd02_x_ind] = 1.3531+1+1+0.5+1
initial[bd02_y_ind] = 0.0

initial[bd03_x_ind] = 1.3531+1
initial[bd03_y_ind] = 1.56+1.56+0.78+1

func(initial)
jacoby(initial)

mass = zeros(number_of_dofs(sys), number_of_dofs(sys));
for i in 1:last_body_dof(sys)
    mass[i, i] = 1
end
time_span = 0:0.005:10
sol1 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol1, initial, mass, func, jacoby, step(time_span))

animate(sys, sol1, time_span, "out/real_5bar_dyn.mp4"; framerate = 60, limits = (-7,7, -7, 7))

sol2 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))

static_solver!( sol2 , initial, func, jacoby)

animate(sys, sol2, time_span, "out/real_5bar_stat.mp4"; framerate = 60, limits = (-7,7, -7, 7))