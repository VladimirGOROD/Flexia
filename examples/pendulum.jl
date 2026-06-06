using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff

using StaticArrays

const g = 9.81

bd1 = Body2D(10, 500)
bd2 = Body2D(10, 500)
bd3 = Body2D(10, 500)
bd4 = Body2D(10, 500)
bd5 = Body2D(10, 500)
bd6 = Body2D(10, 500)
# bd2 = Body2D(10, 1)

bd1.forces[2] = (x, t) -> -bd1.mass * g
bd2.forces[2] = (x, t) -> -bd2.mass * g
bd3.forces[2] = (x, t) -> -bd3.mass * g
bd4.forces[2] = (x, t) -> -bd4.mass * g
bd5.forces[2] = (x, t) -> -bd5.mass * g


jnt1 = FixedJoint(bd1)

jnt2 = HingeJoint(bd1, bd2)
jnt3 = HingeJoint(bd2, bd3)
jnt4 = HingeJoint(bd3, bd4)
jnt5 = HingeJoint(bd4, bd5)
jnt6 = HingeJoint(bd5, bd6)

jnt7 = FixedJoint(bd6)

tcp1 = TorsionalSpring(jnt3, 1100.,0.0, 100.)
tcp2 = TorsionalSpring(jnt4, 1100.,0.0, 100.)
tcp3 = TorsionalSpring(jnt5, 1100.,0.0, 100.)

# traj = TrajectoryJoint(bd3, circular_trajectory([0.4, 1.1], 1.1, 1), 0.0, 10.0)

set_position_on_second_body!(jnt2, SA[-1., 0])

set_position_on_first_body!(jnt3, SA[1., 0])
set_position_on_second_body!(jnt3, SA[-1., 0])
set_position_on_first_body!(jnt4, SA[1., 0])
set_position_on_second_body!(jnt4, SA[-1., 0])
set_position_on_first_body!(jnt5, SA[1., 0])
set_position_on_second_body!(jnt5, SA[-1., 0])

set_position_on_first_body!(jnt6, SA[1., 0])

jnt7.pos = SA[2., 0.0]

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)
add!(sys, bd3)
add!(sys, bd4)
add!(sys, bd5)
add!(sys, bd6)

add!(sys, jnt1)
add!(sys, jnt2)
add!(sys, jnt3)
add!(sys, jnt4)
add!(sys, jnt5)
add!(sys, jnt6)
add!(sys, jnt7)

add!(sys, tcp1)
add!(sys, tcp2)
add!(sys, tcp3)

# add!(sys, traj)
#mb1 = first_marker("mb1", bd1)

if (!assemble!(sys))
    println("Assembling failed!")
end
# end # module Flexia
func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)


bd1_x_ind, bd1_y_ind, _ = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)
bd3_x_ind, bd3_y_ind, bd3_t_ind = get_body_position_dofs(sys, bd3)
bd4_x_ind, bd4_y_ind, bd4_t_ind = get_body_position_dofs(sys, bd4)
bd5_x_ind, bd5_y_ind, bd5_t_ind = get_body_position_dofs(sys, bd5)
bd6_x_ind, bd6_y_ind, _ = get_body_position_dofs(sys, bd6)

bd2_Vx_ind, bd2_Vy_ind, bd2_Vt_ind = get_body_velocity_dofs(sys, bd2)
bd5_Vx_ind, bd5_Vy_ind, bd5_Vt_ind = get_body_velocity_dofs(sys, bd5)

bd3_Vx_ind, bd3_Vy_ind, bd3_Vt_ind = get_body_velocity_dofs(sys, bd3)


initial = zeros(number_of_dofs(sys))
# first configuretion
# initial[bd2_x_ind] = -0.171
# initial[bd2_y_ind] = 0.4698
# initial[bd2_t_ind] = 110*pi/180

# initial[bd3_x_ind] = 0.079
# initial[bd3_y_ind] = 1.2094
# initial[bd3_t_ind] = 32.39*pi/180

# initial[bd4_x_ind] = 0.921
# initial[bd4_y_ind] = 1.2094
# initial[bd4_t_ind] = -32.39*pi/180

# initial[bd5_x_ind] = 1.171
# initial[bd5_y_ind] = 0.4698
# initial[bd5_t_ind] = -110*pi/180

# initial[bd6_x_ind] = 1

 
 initial[bd2_Vt_ind] = -10
 initial[bd5_Vt_ind] = -10

#  zero config
# initial[bd2_x_ind] = 0.5
# initial[bd3_x_ind] = 1
# initial[bd3_y_ind] = 0.5
# initial[bd3_t_ind] = pi/2
# initial[bd4_x_ind] = 1
# initial[bd4_y_ind] = 0.5
# initial[bd4_t_ind] = -pi/2
# initial[bd5_x_ind] = 1.5
# initial[bd6_x_ind] = 2.0

#  second config

initial[bd2_y_ind] = 1
initial[bd2_t_ind] = pi/2

initial[bd3_x_ind] = 0.5
initial[bd3_y_ind] = 2 + sin(pi/3)
initial[bd3_t_ind] = pi/3

initial[bd4_x_ind] = 1.5
initial[bd4_y_ind] = 3
initial[bd4_t_ind] = - pi/3

initial[bd5_x_ind] = 2
initial[bd5_y_ind] = 1
initial[bd5_t_ind] = -pi/2

initial[bd6_x_ind] = 1

func(initial)
jacoby(initial)

mass = get_mass_matrix(sys)
time_span = 0:0.01:10
sol = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol, initial, mass, func, jacoby, step(time_span))

animate(sys, sol, time_span, "out/pendulum.mp4"; framerate = 30, limits = (-5,5, -5, 5))
