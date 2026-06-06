using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff

using StaticArrays
using LinearAlgebra

const g = 9.81

const m1 = 10
const m2 = 20
const m3 = 30

bd1 = Body2D(m1, 100)
bd2 = Body2D(m2, 100)
bd3 = Body2D(m3, 100)

# bd2 = Body2D(10, 1)

bd1.forces[2] = (x, t) -> -m1 * bd1.mass * g
bd2.forces[2] = (x, t) -> -m2 * bd2.mass * g
bd3.forces[2] = (x, t) -> -m3 * bd3.mass * g

jnt1 = FixedJoint(bd1)
jnt2 = HingeJoint(bd1, bd2)
jnt3 = HingeJoint(bd2, bd3)

tcp1 = TorsionalSpring(jnt2, 100000.,0.0, 0.)
tcp2 = TorsionalSpring(jnt3, 100000.,0.0, 0.)

set_position_on_first_body!(jnt2, SA[1.,0])
set_position_on_second_body!(jnt2, SA[-1., 0])

set_position_on_first_body!(jnt3, SA[1., 0])
set_position_on_second_body!(jnt3, SA[-1., 0])


sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)
add!(sys, bd3)


add!(sys, jnt1)
add!(sys, jnt2)
add!(sys, jnt3)

add!(sys, tcp1)
add!(sys, tcp2)

if (!assemble!(sys))
    println("Assembling failed!")
end
# end # module Flexia
func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

bd1_x_ind, bd1_y_ind, _ = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)
bd3_x_ind, bd3_y_ind, bd3_t_ind = get_body_position_dofs(sys, bd3)

bd3_Vx_ind, bd3_Vy_ind, bd3_Vt_ind = get_body_velocity_dofs(sys, bd3)

initial = zeros(number_of_dofs(sys))

initial[bd2_x_ind] = 2.
initial[bd2_y_ind] = 0.
initial[bd3_x_ind] = 3.75
initial[bd3_y_ind] = 0.1
initial[bd3_t_ind] = asin(0.1 / 1)

# initial[bd3_Vt_ind] = 1
func(initial)
jacoby(initial)

mass = get_mass_matrix(sys)

time_start = 0
time_end = 1
time_step = 100
time_span = range(time_start, time_end, time_step)
sol = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
# cros!(sol, initial, mass, func, jacoby, step(time_span))
static_solver!( sol , initial, func, jacoby)
animate(sys, sol, time_span, "out/stat_bar13.mp4"; framerate = floor(Int64, 1.0/step(time_span)), limits = (-5,5, -5, 5))
#
R1 = get_lms(sys, jnt1)
R2 = get_lms(sys, jnt2)
R3 = get_lms(sys, jnt3)

M1 = get_spring_moment(tcp1, 0., sol[bd2_t_ind, time_step])
M2 = get_spring_moment(tcp2, sol[bd2_t_ind, time_step], sol[bd3_t_ind, time_step])

lms11 = sol[R1[1], time_step]
lms12 = sol[R1[2], time_step]
lms13 = sol[R1[3], time_step]

lms21 = sol[R2[1], time_step]
lms22 = sol[R2[2], time_step]

lms31 = sol[R3[1], time_step]
lms32 = sol[R3[2], time_step]

# углы
θ21 = sol[bd2_t_ind, time_step]
θ31 = sol[bd3_t_ind, time_step]

θ2 = rad2deg(sol[bd2_t_ind, time_step])
θ3 = rad2deg(sol[bd3_t_ind, time_step])

println("углы пружин в градусах $θ2, $θ3")
println("углы пружин в радианах $θ21, $θ31")

println("Статический расчёт.")
println("Реакция шарнира 1 по оси Х,Y и момент заделки = [ $lms11 , $lms12, $lms13]")
println("Реакция шарнира 2 по оси Х,Y =  [$lms21, $lms22]")
println("Реакция шарнира 3 по оси Х,Y = [$lms31, $lms32]")

println("Реакция пружины 1 момент = $M1")
println("Реакция пружины 2 момент = $M2")

function write_matrix_formatted(filename::String, matrix::AbstractMatrix)
    open(filename, "w") do file
        # Записываем открывающую скобку
        write(file, "[")
        
        nrows, ncols = size(matrix)
        
        for i in 1:nrows
            # Переходим на новую строку для всех строк кроме первой
            if i > 1
                write(file, " ")
            end
            
            # Записываем элементы строки
            for j in 1:ncols
                write(file, string(matrix[i, j]))
                if j < ncols
                    write(file, " ")
                end
            end
            
            # Записываем разделитель строк или закрывающую скобку
            if i < nrows
                write(file, "\n")
            else
                write(file, "]")
            end
        end
    end
end

initial2 = copy(sol[:,time_step])

initial2[bd3_Vt_ind] = 1

sol2 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))

cros!(sol2, initial2, mass, func, jacoby, step(time_span))

animate(sys, sol2, time_span, "out/dyn_bar13.mp4"; framerate = 30, limits = (-5,5, -5, 5))

lms111 = sol2[R1[1], time_step]
lms122 = sol2[R1[2], time_step]
lms133 = sol2[R1[3], time_step]

lms211 = sol2[R2[1], time_step]
lms222 = sol2[R2[2], time_step]

lms311 = sol2[R3[1], time_step]
lms322 = sol2[R3[2], time_step]

M11 = get_spring_moment(tcp1, 0., sol2[bd2_t_ind, time_step])
M22 = get_spring_moment(tcp2, sol2[bd2_t_ind, time_step], sol2[bd3_t_ind, time_step])

println("\n")

println("Динамический расчёт. Данные в конечный момент времени заданного интервала")
println("Реакция шарнира 1 по оси Х,Y и момент заделки = [ $lms111 , $lms122, $lms133]")
println("Реакция шарнира 2 по оси Х,Y =  [$lms211, $lms222]")
println("Реакция шарнира 3 по оси Х,Y = [$lms311, $lms322]")

println("Реакция пружины 1 момент = $M11")
println("Реакция пружины 2 момент = $M22")

write_matrix_formatted("out/solution.txt", sol)

write_matrix_formatted("out/sol_dynamic.txt", sol2)


f = Figure()

react_x1 = Float64[]
react_x2 = Float64[]
react_x3 = Float64[]

react_y1 = Float64[]
react_y2 = Float64[]
react_y3 = Float64[]

react_M1 = Float64[]

react_TCP1 = Float64[]
react_TCP2 = Float64[]

for i in 1:time_step
    push!(react_x1, sol2[R1[1], i])
    push!(react_x2, sol2[R2[1], i])
    push!(react_x3, sol2[R3[1], i])

    push!(react_y1, sol2[R1[2], i])
    push!(react_y2, sol2[R2[2], i])
    push!(react_y3, sol2[R3[2], i])

    push!(react_M1, sol2[R1[3], i])

    push!(react_TCP1, get_spring_moment(tcp1, 0., sol2[bd2_t_ind, i]))
    push!(react_TCP2, get_spring_moment(tcp2, sol2[bd2_t_ind, i], sol2[bd3_t_ind, i]))
end
react_x1
react_x2
react_x3

react_y1
react_y2

react_M1
time = 1:1:time_step
ax2 = Axis(f[1,3], title="График реакций по оси X",xlabel="Время, с.",ylabel="Реакции, Н.")

l1 = lines!(ax2, time, react_x1, linestyle = :dot)
l2 = lines!(ax2, time, react_x2, linestyle = :dash)
l3 = lines!(ax2, time, react_x3, linestyle = :dashdot)

Legend(f[1 , 4], [l1,l2,l3], ["R1(Fix)", "R2(Jnt1)", "R3(Jnt2)"], framevisible = false, halign = :right, valign = :top)

ax1 = Axis(f[1 , 1], title="График реакций по оси Y",xlabel="Время, с.",ylabel="Реакции, Н.")

l4 = lines!(ax1, time, react_y1, linestyle = :dot)
l5 = lines!(ax1, time, react_y2, linestyle = :dash)
l6 = lines!(ax1, time, react_y3)

Legend(f[1 , 2], [l4,l5,l6], ["R1(Fix)", "R2(Jnt1)", "R3(Jnt2)"], framevisible = false, halign = :left, valign = :top)

ax3 = Axis(f[2,1], title="График моментов пружин и заделки",xlabel="Время, с.",ylabel="Момент, Нxм.")

l7 = lines!(ax3, time, react_M1)
l8 = lines!(ax3, time, react_TCP1, linestyle = :dot)
l9 = lines!(ax3, time, react_TCP2, linestyle = :dash)

Legend(f[2 , 2], [l7,l8,l9], ["M1(Fix)", "TCP1", "TCP2"], framevisible = false, halign = :left, valign = :top)

save("out/Reactions.png", f)
