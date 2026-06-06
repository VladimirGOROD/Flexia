using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using StaticArrays

const g = 9.81

m = 1; M = 5; L = 1.; d = 1.;
ground = Body2D(1e6, 1e6; length=0.01)  # Large mass/inertia to simulate fixed ground
crank = Body2D(m, 0.1; length=L)
slider = Body2D(M, 0.1; length=0.25)

crank.forces[2] = (x, t) -> -crank.mass * g
slider.forces[2] = (x, t) -> -slider.mass * g

sys = MBSystem2D()
add!(sys, ground)
add!(sys, crank)
add!(sys, slider)

ground_joint = FixedJoint(ground)
setposition!(ground_joint, SA[0.0, 0.0])
setrotation!(ground_joint, 0.0)

slider_joint = SliderJoint(ground, slider)
set_position_on_first_body!(slider_joint, SA[0.0, 0.0])
set_position_on_second_body!(slider_joint, SA[0.0, 0.0])
set_direction_on_first_body!(slider_joint, SA[1.0, 0.0])
set_direction_on_second_body!(slider_joint, SA[1.0, 0.0])

hinge_joint = HingeJoint(slider, crank)
set_position_on_first_body!(hinge_joint, SA[0.0, 0.0])
set_position_on_second_body!(hinge_joint, SA[-crank.length/2, 0.0])

slider_actuator = PositionLinearActuator2D(slider_joint, 0.0)

spring_damper = TorsionalSpring(hinge_joint; stiffness=0.0, damping=0.1, vis_r=0.1)

# Add joints to system
add!(sys, ground_joint)
add!(sys, slider_joint)
add!(sys, hinge_joint)

add!(sys, slider_actuator)
add!(sys, spring_damper)

initial_state = zeros(number_of_dofs(sys))

crank_θ = π/2
crank_x = -1 + crank.length/2 * cos(crank_θ)
crank_y = crank.length/2 * sin(crank_θ)
set_initial_position!(initial_state, sys, crank, SA_F64[crank_x, crank_y, crank_θ])

slider_x = -1
set_initial_position!(initial_state, sys, slider, SA_F64[slider_x, 0, 0])
assemble!(sys)

# set_initial_velocity!(initial_state, sys, crank, SA_F64[0, -5, 0])
#

measured = [
    get_body_position_dofs(sys, slider)[1], 
    get_body_velocity_dofs(sys, slider)[1], 
    get_body_position_dofs(sys, crank)[3], 
    get_body_velocity_dofs(sys, crank)[3]
    ]

C = zeros(length(measured), number_of_dofs(sys)-1)
for i in eachindex(measured)
    C[i, measured[i]] = 1
end

controlled = get_lms(sys, slider_actuator)
B = zeros(number_of_dofs(sys)-1, length(controlled))
for i in eachindex(controlled)
    B[controlled[i], i] = 1
end

A = sys.jacobian(initial_state)[1:end-1, 1:end-1]


# K = lqr(A, B, C, 0.001)

K = 0.01 * SA_F64[1 1 -1 -1]
ŷ = SA_F64[1, 0, π/2, 0]
K * ŷ
sys.prestep = (state) -> begin
    y = state[measured];
    t = state[end]
    # err = ŷ - y
    # u = K * (ŷ - y)
    # u = y[1] + 0.005 * err[1]
    settarget!(slider_actuator,  slider_x + sin(t - π/2) + sin(π/2))
end

crank_θ = π/2 + 0.0001

#
#
time_span = 0:0.0025:10
ratio = 1
sol = simulate(sys, initial_state, time_span)

animate(sys, sol, ratio*time_span, "out/inverted_pendulum.mp4"; framerate=floor(Int64, 1.0/(ratio * step(time_span))), limits=(-2, 2, -2, 2))

##
for_spy(mat) = reverse(mat', dims=2)
spy(for_spy(sys.jacobian(initial_state)), axis=(;aspect=DataAspect()))