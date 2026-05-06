using GLMakie

function legth_param(l_1::Float64,l_2::Float64, l_3::Float64, l_4::Float64, l_5::Float64)
    a1 = l_1
    a2 = l_3
    b1 = l_2
    b2 = l_4
    w = l_5
    return a1, a2, b1, b2, w
end
a1, a2, b1, b2, w = legth_param(1.0, 1.0, 1.0 , 1.0, 1.0)

function inverse_kinematics(x_P::Float64, y_P::Float64)

    c1 = sqrt(x_P^2 + y_P^2)
    c2 = sqrt(c1^2 - 2 * x_P * w + w^2)
    
    α1 = acosd(x_P / c1)
    α2 = acosd( (-x_P + w) / c2)

    β1 = acosd((b1^2 - a1^2 - c1^2 ) / ( -2 * a1 * c1))
    β2 = acosd((b2^2 - a2^2 - c2^2 ) / (-2 * a2 * c2))

    θ1 = α1 + β1
    θ2 = α2 + β2

    return θ1, θ2

end

function forward_kinematics(θ1::Float64, θ2::Float64)
    x_R1 = a1 * cosd(θ1)
    y_R1 = a1 * sind(θ1)

    x_R2 = a2 * cosd(θ2) + w
    y_R2 = a2 * sind(θ2)
# лямбды при вычитании двух окружных уравнений друг из друга
    λ1 = (y_R1 - y_R2) / (x_R2 - x_R1)
    λ2 = (b1^2 - b2^2 - a1^2 + x_R2^2 + y_R2^2) / (2 * (x_R2 - x_R1))

# лямюды при подстановки линейной комбинации x_P(y_P) и решение квадратного уравнения
    λ3 = 1 + λ1^2 # a - в квадратном ур-е
    λ4 = 2 * (λ1 * λ2 - λ1 * x_R1 - y_R1) # b - в квадратном ур-е
    λ5 = a1^2 - b1^2 + λ2^2 - 2 * λ2 * x_R1 # c - в квадратном ур-е
   
    x_P1 = 0
    y_P1 = 0

    x_P2 = 0
    y_P2 = 0

    if (λ4^2 - 4 * λ3 * λ5) >=0   #условие избигания критических положений, где дискриминант равен sqrt(<0)
# находим по дискриминанту два решения
        y_P1 = (-λ4 + sqrt(λ4^2 - 4 * λ3 * λ5)) / (2 * λ3)
        y_P2 = (-λ4 - sqrt(λ4^2 - 4 * λ3 * λ5)) / (2 * λ3)

        x_P1 = y_P1 * λ1 + λ2
        x_P2 = y_P2 * λ1 + λ2
    end
    return [x_P1, y_P1, x_P2, y_P2, x_R1, y_R1, x_R2, y_R2]
end

x_traj = Float64[]
y_traj = Float64[]

x_traj = [ 0.3, 0.4 , 0.5, 0.6]

y_traj = [1.866, 1.866, 1.866, 1.866]
f = Figure()
ax = Axis(f[1,1], xlabel="X координата", ylabel="Y координата", aspect = DataAspect())
lines!(ax, x_traj, y_traj, color = :black)
cs = rand(100)
function draw_mechanism_on_traj(f, ax, x_traj, y_traj)
    x_mech = Float64[]
    y_mech = Float64[]
    for x in x_traj
        for y in y_traj
            τ1 , τ2 = inverse_kinematics(x, y)
            p_x_mech = [0. , a1 * cosd(τ1), x, a2 * cosd(180 - τ2) + w, w]
            p_y_mech = [0. , a1 * sind(τ1), y, a2 * sind(180 - τ2) , 0.]

            append!(x_mech, p_x_mech)
            append!(y_mech, p_y_mech)
        end
    end
    lines!(ax, x_mech, y_mech, color = cs)    
end
draw_mechanism_on_traj(f, ax, x_traj, y_traj)
save("mechanism_on_traj.png", f)

n_points = 500  # Количество точек по каждому углу
θ_min = 0    # Минимальное значение угла
θ_max = 360      # Максимальное значение угла


θ1_values = range(θ_min, θ_max, length=n_points)
θ2_values = range(θ_min, θ_max, length=n_points)


all_points_x = Float64[]
all_points_y = Float64[]

# закрашивание рабочей области построением механизма во всех положения
function draw_workspace!(f,ax)
    for θ1 in θ1_values
        for θ2 in θ2_values
            result2 = forward_kinematics(θ1, θ2)

            all_x = [0, result2[5], result2[1], result2[7], result2[3], result2[5], result2[3], result2[7], w]
            all_y = [0, result2[6], result2[2], result2[8], result2[4], result2[6], result2[4], result2[8], 0]

            append!(all_points_x, all_x)
            append!(all_points_y, all_y)
        end
    end
    lines!(ax, all_points_x, all_points_y, color=:orange)
end
ax = Axis(f[1, 2], title="Рабочая область механизма",xlabel="X координата",ylabel="Y координата",aspect=DataAspect())
draw_workspace!(fig,ax)
save("mechamism's_workspace.png", f)