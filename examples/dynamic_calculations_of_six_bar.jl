using Pkg; Pkg.activate("./examples")
using Flexia
using ForwardDiff
using GLMakie
using StaticArrays
# using CSV
# using DataFrames
const g = 9.81

CURRENT_TIME = 0.

# K₁ = 0.03
# K₂ = 0.0115
K_hsp = 500.
K₁ = 0.5
K₂ = K₁
K₃ = K₂
K₄ = K₂
K₅ = K₂
K₆ = K₁


bd2_bd6_offset_x = 0.312
actuators_offset_y = 0.1084
const angle_5 = deg2rad(-110)

# Массы в кг
const m1 = 24e-3
const m2 = 29e-3
const m4 = 51e-3

const m5 = 52e-3
# Моменты инерции в кг·м²
# 104726.798394 г·мм² = 104726.798394e-9 кг·м²
I2 = 104726.798394e-9
# 28187.202983 г·мм² = 28187.202983e-9 кг·м²
I4 = 36119.84e-9
F_container = 1
omega = 0.5
# Длины в метрах: было 1.5 дм = 0.15 м, 0.75 дм = 0.075 м, 0.5 дм = 0.05 м
bd1 = Body2D(m5, 0.1e-30, length = 0.15)
bd2 = Body2D(m2+m1, I2, length = 0.15)
bd3 = Body2D(m2+m1, I2, length = 0.15)

bd4 = Body2D(m4+m1, I4, length = 0.1)
# bd4.forces[1] = (x) -> -(m4+m1) * F_container
# bd4.forces[1] = (x, t) -> F_container * sin(omega * t) +rand()*0.01
bd5 = Body2D(m2+m1, I2, length = 0.15)
bd6 = Body2D(m2+m1, I2, length = 0.15)
bd7 = Body2D(m5, 0.1e-30, length = 0.15)

slider_rail_1 = Body2D(m2+m1, I2, length = 0.068)
slider_rail_2 = Body2D(m2+m1, I2, length = 0.068)

slider_1 = Body2D(m4, I4, length = 0.0255)
slider_2 = Body2D(m4, I4, length = 0.0255)

jnt1 = FixedJoint(bd1)
jnt8 = FixedJoint(bd7)

jnt2 = HingeJoint(bd1, bd2)
jnt3 = HingeJoint(bd2, bd3)

jnt4 = HingeJoint(bd3, bd4)
jnt5 = HingeJoint(bd4, bd5)

jnt6 = HingeJoint(bd5, bd6)
jnt7 = HingeJoint(bd6, bd7)

# Left side with actuator

hinge_x_offset = 0.0
ground_rail_hinge1 = HingeJoint(bd1, slider_rail_1)
set_position_on_first_body!(ground_rail_hinge1, SA[hinge_x_offset, actuators_offset_y])
set_position_on_second_body!(ground_rail_hinge1, SA[0., 0.])


direcrion_SL = SA[-1., 0.]

x_offset = bd1.length/2 - slider_1.length/2

slider_ground_slider_1 = SliderJoint(slider_rail_1, slider_1)
set_position_on_first_body!(slider_ground_slider_1, SA[x_offset, 0.])
set_position_on_second_body!(slider_ground_slider_1, SA[0., 0.])
set_direction_on_first_body!(slider_ground_slider_1, direcrion_SL)
set_direction_on_second_body!(slider_ground_slider_1, direcrion_SL)

hor_bd2_hinge = HingeJoint(slider_1, bd2)
set_position_on_first_body!(hor_bd2_hinge, SA[slider_1.length/2, 0])
set_position_on_second_body!(hor_bd2_hinge, SA[-(bd2.length/2-actuators_offset_y), 0])


# right side with actuator

ground_rail_hinge2 = HingeJoint(bd7, slider_rail_2)
set_position_on_first_body!(ground_rail_hinge2, SA[0., actuators_offset_y])
set_position_on_second_body!(ground_rail_hinge2, SA[0., 0.])

hor_bd6_hinge = HingeJoint(bd6, slider_2)
set_position_on_first_body!(hor_bd6_hinge, SA[-(actuators_offset_y-bd6.length/2), 0])
set_position_on_second_body!(hor_bd6_hinge, SA[-(slider_2.length/2), 0])

direcrion_SL_rev = SA[-1., 0.]
slider_ground_slider_2 = SliderJoint(slider_rail_2, slider_2)
set_position_on_first_body!(slider_ground_slider_2, SA[-x_offset, 0.])
set_position_on_second_body!(slider_ground_slider_2, SA[0., 0.])
set_direction_on_first_body!(slider_ground_slider_2, direcrion_SL_rev)
set_direction_on_second_body!(slider_ground_slider_2, direcrion_SL_rev)

tcp1 = TorsionalSpring(jnt2, K₁, deg2rad(-90), 0.1, 0.03)
tcp2 = TorsionalSpring(jnt3, K₂, deg2rad(45), 0.1, 0.03)

tcp3 = TorsionalSpring(jnt4, K₃, deg2rad(45), 0.1, 0.03)
tcp4 = TorsionalSpring(jnt5, K₄, deg2rad(45), 0.1, 0.03)

tcp5 = TorsionalSpring(jnt6, K₅, deg2rad(45), 0.1, 0.03)
tcp6 = TorsionalSpring(jnt7, K₆, deg2rad(-90), 0.1, 0.03)

tcp_ground_rail_1 = TorsionalSpring(ground_rail_hinge1, K₁, deg2rad(0), 0.1, 0.03)
tcp_ground_rail_2 = TorsionalSpring(ground_rail_hinge2, K₁, deg2rad(0), 0.1, 0.03)

tcp_hor_bd2 = TorsionalSpring(hor_bd2_hinge, K₁, deg2rad(-90), 0.1, 0.03)
tcp_hor_bd6 = TorsionalSpring(hor_bd6_hinge, K₁, deg2rad(-90), 0.1, 0.03)

hsp1 = LinearSpring(slider_ground_slider_1;
                    stiffness = K_hsp,
                    damping = 0.01, vis_r = 0.1, vis_N = 6)
hsp2 = LinearSpring(slider_ground_slider_2;
                    stiffness = K_hsp, 
                    damping = 0.01, vis_r = 0.1, vis_N = 6)

# Позиции присоединений: всё было в дм, теперь в м (делим на 10)
set_position_on_first_body!(jnt2, SA[bd1.length/2, 0.])
set_position_on_second_body!(jnt2, SA[-bd2.length/2, 0.])

set_position_on_first_body!(jnt3, SA[bd2.length / 2, 0.])
set_position_on_second_body!(jnt3, SA[-bd3.length /2, 0.])

set_position_on_first_body!(jnt4, SA[bd3.length/2, 0])
set_position_on_second_body!(jnt4, SA[-bd4.length/2, 0])

set_position_on_first_body!(jnt5, SA[bd4.length/2, 0])
set_position_on_second_body!(jnt5, SA[-bd5.length/2, 0])

set_position_on_first_body!(jnt6, SA[bd5.length/2, 0])
set_position_on_second_body!(jnt6, SA[-bd6.length/2, 0])

set_position_on_first_body!(jnt7, SA[bd6.length/2, 0])
set_position_on_second_body!(jnt7, SA[-bd7.length/2, 0])


# Позиция фиксации тела 7: 6.12 дм = 0.612 м
jnt8.pos = SA[bd1.length/2 + bd2_bd6_offset_x + bd7.length/2, 0.]

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)
add!(sys, bd3)
add!(sys, bd4)
add!(sys, bd5)
add!(sys, bd6)
add!(sys, bd7)

add!(sys, slider_rail_1)
add!(sys, slider_1)

add!(sys, slider_rail_2)
add!(sys, slider_2)

add!(sys, jnt1)

add!(sys, jnt2)
add!(sys, jnt3)
add!(sys, jnt4)
add!(sys, jnt5)
add!(sys, jnt6)
add!(sys, jnt7)

add!(sys, jnt8)

add!(sys, ground_rail_hinge1)
add!(sys, ground_rail_hinge2)

add!(sys, slider_ground_slider_1)
add!(sys,slider_ground_slider_2)

add!(sys,hor_bd2_hinge)
add!(sys,hor_bd6_hinge)

add!(sys, tcp1)
add!(sys, tcp2)
add!(sys, tcp3)
add!(sys, tcp4)
add!(sys, tcp5)
add!(sys, tcp6)

add!(sys,tcp_ground_rail_1)
add!(sys,tcp_ground_rail_2)

add!(sys,tcp_hor_bd2)
add!(sys,tcp_hor_bd6)

add!(sys,hsp1)
add!(sys,hsp2)

if (!assemble!(sys))
    println("Assembling failed!")
end

func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

bd1_x_ind, bd1_y_ind, bd1_t_ind = get_body_position_dofs(sys, bd1)

bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)

bd3_x_ind, bd3_y_ind, bd3_t_ind = get_body_position_dofs(sys, bd3)

bd4_x_ind, bd4_y_ind, bd4_t_ind = get_body_position_dofs(sys, bd4)
bd4_vx_ind, bd4_vy_ind, bd4_vt_ind = get_body_velocity_dofs(sys, bd4)

bd5_x_ind, bd5_y_ind, bd5_t_ind = get_body_position_dofs(sys, bd5)

bd6_x_ind, bd6_y_ind, bd6_t_ind = get_body_position_dofs(sys, bd6)

bd7_x_ind, bd7_y_ind, bd7_t_ind = get_body_position_dofs(sys, bd7)

bd4_Vx_ind, bd4_Vy_ind, bd4_Vt_ind = get_body_velocity_dofs(sys, bd4)


initial = zeros(number_of_dofs(sys))


# Все начальные координаты переведены из дм в м (делением на 10)
initial[bd1_x_ind] = 0.

initial[bd2_x_ind] = bd1.length/2
initial[bd2_y_ind] = bd2.length/2
initial[bd2_t_ind] = deg2rad(90)

initial[bd3_x_ind] = bd1.length/2 + bd3.length*cos(deg2rad(45))/2
initial[bd3_y_ind] = bd2.length + bd3.length*sin(deg2rad(45))/2
initial[bd3_t_ind] = deg2rad(45)

initial[bd4_x_ind] = bd1.length/2 + bd3.length*cos(deg2rad(45)) + bd4.length/2
initial[bd4_y_ind] = bd2.length + bd3.length*sin(deg2rad(45))
initial[bd4_t_ind] = deg2rad(0)

initial[bd5_x_ind] = bd1.length/2 + bd2_bd6_offset_x - bd5.length*cos(deg2rad(45))/2
initial[bd5_y_ind] = bd6.length + bd5.length*sin(deg2rad(45))/2
initial[bd5_t_ind] = deg2rad(-45)

initial[bd6_x_ind] = bd1.length/2 + bd2_bd6_offset_x
initial[bd6_y_ind] = bd6.length/2
initial[bd6_t_ind] = deg2rad(-90)

slider_rail_1_X, slider_rail_1_Y, _ =  get_body_position_dofs(sys, slider_rail_1)
initial[slider_rail_1_X] = 0.
initial[slider_rail_1_Y] = actuators_offset_y

slider_rail_2_X, slider_rail_2_Y, _ =  get_body_position_dofs(sys, slider_rail_2)
initial[slider_rail_2_X] = bd1.length/2 + bd2_bd6_offset_x + bd7.length/2
initial[slider_rail_2_Y] = actuators_offset_y



slider_1_X, slider_1_Y, _ =  get_body_position_dofs(sys, slider_1)
initial[slider_1_X] = bd1.length/2 - slider_1.length/2
initial[slider_1_Y] = actuators_offset_y

slider_2_X, slider_2_Y, _ =  get_body_position_dofs(sys, slider_2)
initial[slider_2_X] = bd1.length/2 + bd2_bd6_offset_x + slider_2.length/2
initial[slider_2_Y] = actuators_offset_y

# Начальное время (последний элемент)
initial[end] = 0.0

initial[bd7_x_ind] = bd1.length/2 + bd2_bd6_offset_x + bd7.length/2


f = (t)-> func([initial[1:end-1]; t])[bd4_vx_ind]



jacoby(initial)
#

mass = get_mass_matrix(sys)
t_end = 10*π / omega
time_span = LinRange(0,t_end, 501)

n_steps = length(time_span)

sol1 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol1, initial, mass, func, jacoby, step(time_span))
# lines(sol1[end, :])
# Лимиты графика тоже в метрах
animate(sys, sol1, time_span, "Flexia/out/new_6bar.mp4"; framerate = 30, limits = (-0.1, 0.6, -0.1, 0.5))

# Flexia.draw_static(sys, initial; limits = (-0.1, 0.6, -0.1, 0.5))

##
# fig = Figure()
# ax = Axis(fig[1,1], aspect = DataAspect())
# for body in Flexia.bodies(sys)
#     bar = Vector{Point2f}(undef, 2)
#     bar2 = Vector{Point2f}(undef, 2)
#     bar .= get_boundary_points(sys, body, initial)
#     bar2 .= get_boundary_points(sys, body, sol1[:,2])
#     lines!(ax, bar)
#     lines!(ax, bar2, linestyle=:dash)
# end
# fig
sys
