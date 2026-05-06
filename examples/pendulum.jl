using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff

using StaticArrays

const g = 9.81

bd1 = Body2D(10, 1)
bd2 = Body2D(10, 1)
bd3 = Body2D(10, 1)
bd4 = Body2D(10, 1)
bd5 = Body2D(10, 1)
# bd2 = Body2D(10, 1)

bd1.forces[2] = (x) -> -bd1.mass * g
bd2.forces[2] = (x) -> -bd2.mass * g
bd3.forces[2] = (x) -> -bd3.mass * g
bd4.forces[2] = (x) -> -bd4.mass * g
bd5.forces[2] = (x) -> -bd5.mass * g


jnt1 = FixedJoint(bd1)
jnt2 = HingeJoint(bd1, bd2)
jnt3 = HingeJoint(bd2, bd3)
jnt4 = HingeJoint(bd3, bd4)
jnt5 = HingeJoint(bd4, bd5)
set_position_on_second_body!(jnt2, SA[-0.5, 0])

set_position_on_first_body!(jnt3, SA[0.5, 0])
set_position_on_second_body!(jnt3, SA[-0.5, 0])
set_position_on_first_body!(jnt4, SA[0.5, 0])
set_position_on_second_body!(jnt4, SA[-0.5, 0])
set_position_on_first_body!(jnt5, SA[0.5, 0])
set_position_on_second_body!(jnt5, SA[-0.5, 0])

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)
add!(sys, bd3)
add!(sys, bd4)
add!(sys, bd5)

add!(sys, jnt1)
add!(sys, jnt2)
add!(sys, jnt3)
add!(sys, jnt4)
add!(sys, jnt5)

if (!assemble!(sys))
    println("Assembling failed!")
end
# end # module Flexia
func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)


bd1_x_ind, bd1_y_ind, _ = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, _ = get_body_position_dofs(sys, bd2)
bd3_x_ind, bd3_y_ind, _ = get_body_position_dofs(sys, bd3)

initial = zeros(number_of_dofs(sys))
initial[bd2_x_ind] = 0.5
func(initial)
jacoby(initial)

mass = zeros(number_of_dofs(sys), number_of_dofs(sys));
for i in 1:last_body_dof(sys)
    mass[i, i] = 1
end
time_span = 0:0.01:10
sol = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol, initial, mass, func, jacoby, step(time_span))

animate(sys, sol, time_span, "out/pendulum.mp4"; framerate = 30, limits = (-5,5, -5, 5))
