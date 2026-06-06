using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using JSON

using StaticArrays

const g = 9.81

const m1 = 10
const m2 = 10
const m3 = 10
const m4 = 10
const m5 = 10
const m6 = 10
const m7 = 10
const m8 = 10

bd1 = Body2D(m1, 500)
bd2 = Body2D(m2, 500)
bd3 = Body2D(m3, 500)
bd4 = Body2D(m4, 500)
bd5 = Body2D(m5, 500)
bd6 = Body2D(m6, 500)
bd7 = Body2D(m7, 500)
bd8 = Body2D(m8, 500)
# bd2 = Body2D(10, 1)

bd1.forces[2] = (x, t) -> -bd1.mass * g
bd2.forces[2] = (x, t) -> -bd2.mass * g
bd3.forces[2] = (x, t) -> -bd3.mass * g
bd4.forces[2] = (x, t) -> -bd4.mass * g
bd5.forces[2] = (x, t) -> -bd5.mass * g
bd6.forces[2] = (x, t) -> -bd6.mass * g
bd7.forces[2] = (x, t) -> -bd7.mass * g
bd8.forces[2] = (x, t) -> -bd8.mass * g

jnt1 = FixedJoint(bd1)

jnt2 = HingeJoint(bd1,bd2)
jnt3 = HingeJoint(bd2,bd3)
jnt4 = HingeJoint(bd3,bd4)
jnt5 = HingeJoint(bd4,bd5)
jnt6 = HingeJoint(bd5,bd6)
jnt7 = HingeJoint(bd6,bd7)
jnt8 = HingeJoint(bd7,bd8)

jnt9 = FixedJoint(bd8)

tcp2 = TorsionalSpring(jnt3, 100000.,0.0, 0.)
tcp3 = TorsionalSpring(jnt4, 100000.,0.0, 0.)
tcp4 = TorsionalSpring(jnt5, 100000.,0.0, 0.)
tcp5 = TorsionalSpring(jnt6, 100000.,0.0, 0.)
tcp6 = TorsionalSpring(jnt7, 100000.,0.0, 0.)


set_position_on_second_body!(jnt2, SA[-1., 0])

set_position_on_first_body!(jnt3, SA[1., 0])
set_position_on_second_body!(jnt3, SA[-1., 0])

set_position_on_first_body!(jnt4, SA[1., 0])
set_position_on_second_body!(jnt4, SA[-1., 0])

set_position_on_first_body!(jnt5, SA[1., 0])
set_position_on_second_body!(jnt5, SA[-1., 0])

set_position_on_first_body!(jnt6, SA[1., 0])
set_position_on_second_body!(jnt6, SA[-1., 0])

set_position_on_first_body!(jnt7, SA[1., 0])
set_position_on_second_body!(jnt7, SA[-1., 0])

set_position_on_first_body!(jnt8, SA[1., 0])

setposition!(jnt9, SA[7.,7.])
setrotation!(jnt9, pi/2)

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)
add!(sys, bd3)
add!(sys, bd4)
add!(sys, bd5)
add!(sys, bd6)
add!(sys, bd7)
add!(sys, bd8)

add!(sys, jnt1)
add!(sys, jnt2)
add!(sys, jnt3)
add!(sys, jnt4)
add!(sys, jnt5)
add!(sys, jnt6)
add!(sys, jnt7)
add!(sys, jnt8)
add!(sys, jnt9)

add!(sys, tcp2)
add!(sys, tcp3)
add!(sys, tcp4)
add!(sys, tcp5)
add!(sys, tcp6)


if (!assemble!(sys))
    println("Assembling failed!")
end

func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

bd1_x_ind, bd1_y_ind, _ = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)
bd3_x_ind, bd3_y_ind,  = get_body_position_dofs(sys, bd3)
bd4_x_ind, bd4_y_ind, bd4_t_ind = get_body_position_dofs(sys, bd4)
bd5_x_ind, bd5_y_ind, _ = get_body_position_dofs(sys, bd5)
bd6_x_ind, bd6_y_ind, bd6_t_ind = get_body_position_dofs(sys, bd6)
bd7_x_ind, bd7_y_ind, _ = get_body_position_dofs(sys, bd5)
bd8_x_ind, bd8_y_ind, bd8_t_ind = get_body_position_dofs(sys, bd6)

initial = zeros(number_of_dofs(sys))

initial[bd2_x_ind] = 1.
initial[bd2_y_ind] = 1.
initial[bd2_t_ind] = pi/2

initial[bd3_x_ind] = 2.
initial[bd3_y_ind] = 2.

initial[bd4_x_ind] = 3.
initial[bd4_y_ind] = 3.
initial[bd4_t_ind] = pi/2

initial[bd5_x_ind] = 4.
initial[bd5_y_ind] = 4.

initial[bd6_x_ind] = 5.
initial[bd6_y_ind] = 5.
initial[bd6_t_ind] = pi/2

initial[bd7_x_ind] = 6.
initial[bd7_y_ind] = 6.

initial[bd8_x_ind] = 7.
initial[bd8_y_ind] = 7.
initial[bd8_t_ind] = pi/2


func(initial)
jacoby(initial)

mass = get_mass_matrix(sys)

time_span = 0:0.05:200

sol1 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol1, initial, mass, func, jacoby, step(time_span))
animate(sys, sol1, time_span, "out/diag_bar.mp4"; framerate = 60, limits = (-2,13, -2, 13))

sol2 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
static_solver!( sol2 , initial, func, jacoby)
animate(sys, sol2, time_span, "out/diag_bar_stat.mp4"; framerate = 60, limits = (-2, 13, -2, 13))
