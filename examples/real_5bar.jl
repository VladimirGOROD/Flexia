using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff

using StaticArrays

const g = 9.81

# K₁ = 0.03
# K₂ = 0.0115
K₁ = 0.02
K₂ = 0.02
K₃ = K₂
K₄ = K₂
K₅ = K₂
K₆ = K₁


const angle_5 = deg2rad(-90)

# Массы в кг
m1 = 23e-3
m2 = 29e-3
m4 = 15e-3

m5 = 52e-3

# Моменты инерции в кг·м²
# 104726.798394 г·мм² = 104726.798394e-9 кг·м²
I2 = 104726.798394e-9
# 28187.202983 г·мм² = 28187.202983e-9 кг·м²
I4 = 28187.202983e-9


# Длины в метрах: было 1.5 дм = 0.15 м, 0.75 дм = 0.075 м, 0.5 дм = 0.05 м
bd1 = Body2D(m5, 0.1e-30, length = 0.15)
bd2 = Body2D(m2+m1, I2, length = 0.075)
bd3 = Body2D(m2+m1, I2, length = 0.075)

bd4 = Body2D(m4+m1, I4, length = 0.05)

bd5 = Body2D(m2+m1, I2, length = 0.075)
bd6 = Body2D(m2+m1, I2, length = 0.075)
bd7 = Body2D(m5, 0.1e-30, length = 0.15)

jnt1 = FixedJoint(bd1)
jnt8 = FixedJoint(bd7)

jnt2 = HingeJoint(bd1, bd2)
jnt3 = HingeJoint(bd2, bd3)

jnt4 = HingeJoint(bd3, bd4)
jnt5 = HingeJoint(bd4, bd5)

jnt6 = HingeJoint(bd5, bd6)
jnt7 = HingeJoint(bd6, bd7)

tcp1 = TorsionalSpring(bd1, bd2, K₁, deg2rad(-90), 0., 0.015)
tcp2 = TorsionalSpring(bd2, bd3, K₂, deg2rad(45), 0., 0.015)

tcp3 = TorsionalSpring(bd3, bd4, K₃, deg2rad(45), 0., 0.015)
tcp4 = TorsionalSpring(bd4, bd5, K₄, deg2rad(45), 0., 0.015)

tcp5 = TorsionalSpring(bd5, bd6, K₅, deg2rad(45), 0., 0.015)
tcp6 = TorsionalSpring(bd6, bd7, K₆, angle_5, 0., 0.015)

# Позиции присоединений: всё было в дм, теперь в м (делим на 10)
set_position_on_first_body!(jnt2, SA[0.15, 0.])
set_position_on_second_body!(jnt2, SA[-0.075, 0.])

set_position_on_first_body!(jnt3, SA[0.075, 0])
set_position_on_second_body!(jnt3, SA[-0.075, 0])

set_position_on_first_body!(jnt4, SA[0.075, 0])
set_position_on_second_body!(jnt4, SA[-0.05, 0])

set_position_on_first_body!(jnt5, SA[0.05, 0])
set_position_on_second_body!(jnt5, SA[-0.075, 0])

set_position_on_first_body!(jnt6, SA[0.075, 0])
set_position_on_second_body!(jnt6, SA[-0.075, 0])

set_position_on_first_body!(jnt7, SA[0.075, 0])
set_position_on_second_body!(jnt7, SA[-0.15, 0])


# Позиция фиксации тела 7: 6.12 дм = 0.612 м
jnt8.pos = SA[0.612, 0.]

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)
add!(sys, bd3)
add!(sys, bd4)
add!(sys, bd5)
add!(sys, bd6)
add!(sys, bd7)

add!(sys, jnt1)

add!(sys, jnt2)
add!(sys, jnt3)
add!(sys, jnt4)
add!(sys, jnt5)
add!(sys, jnt6)
add!(sys, jnt7)

add!(sys, jnt8)

add!(sys, tcp1)
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

bd1_x_ind, bd1_y_ind, bd1_t_ind = get_body_position_dofs(sys, bd1)

bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)

bd3_x_ind, bd3_y_ind, bd3_t_ind = get_body_position_dofs(sys, bd3)

bd4_x_ind, bd4_y_ind, bd4_t_ind = get_body_position_dofs(sys, bd4)

bd5_x_ind, bd5_y_ind, bd5_t_ind = get_body_position_dofs(sys, bd5)

bd6_x_ind, bd6_y_ind, bd6_t_ind = get_body_position_dofs(sys, bd6)

bd7_x_ind, bd7_y_ind, bd7_t_ind = get_body_position_dofs(sys, bd7)

initial = zeros(number_of_dofs(sys))

# Все начальные координаты переведены из дм в м (делением на 10)
initial[bd1_x_ind] = 0.

initial[bd2_x_ind] = 0.15
initial[bd2_y_ind] = 0.15/2
initial[bd2_t_ind] = deg2rad(90)

initial[bd3_x_ind] = 0.15 + 0.15*cos(deg2rad(45))/2
initial[bd3_y_ind] = 0.15 + 0.15*sin(deg2rad(45))/2
initial[bd3_t_ind] = deg2rad(45)

initial[bd4_x_ind] = 0.15 + 0.05 + 0.15*cos(deg2rad(45))
initial[bd4_y_ind] = 0.15 + 0.15*sin(deg2rad(45))
initial[bd4_t_ind] = deg2rad(0)

initial[bd5_x_ind] = 0.462 - 0.15*cos(deg2rad(46))/2
initial[bd5_y_ind] = 0.15 + 0.15*sin(deg2rad(45))/2
initial[bd5_t_ind] = deg2rad(-45)

initial[bd6_x_ind] = 0.462
initial[bd6_y_ind] = 0.15/2
initial[bd6_t_ind] = deg2rad(-90)

initial[bd7_x_ind] = 0.612

func(initial)
jacoby(initial)

mass = Flexia.get_mass_matrix(sys)
# ##
# time_span = 0:0.005:10
# sol1 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
# cros!(sol1, initial, mass, func, jacoby, step(time_span))

# # Лимиты графика тоже в метрах
# animate(sys, sol1, time_span, "triv_dyn.mp4"; framerate = 60, limits = (-0.1, 0.5, -0.1, 0.5))

