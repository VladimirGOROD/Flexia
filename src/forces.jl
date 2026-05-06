mutable struct TorsionalSpring <: AbstractForce2D
    body1::Body2D
    body2::Body2D
    stiffness::Float64
    rest_angle::Float64
    damping::Float64
    
    function TorsionalSpring(body1::Body2D, body2::Body2D, stiffness::Float64=1.0, rest_angle::Float64=0.0, damping::Float64=0.0)
        return new(body1, body2, stiffness, rest_angle, damping)
    end
end

number_of_dofs(::TorsionalSpring) = 0  # пружина не добавляет лагранжевых множителей

function add!(sys::MBSystem2D, spring::TorsionalSpring)
    push!(sys.forces, spring)
end

function add_force_to_rhs!(rhs, state, sys::MBSystem2D, spring::TorsionalSpring)
    bd1 = spring.body1
    bd2 = spring.body2
    
    # Получаем индексы DOF
    bd1_pos_dofs = get_body_position_dofs(sys, bd1)
    bd1_vel_dofs = get_body_velocity_dofs(sys, bd1)
    bd2_pos_dofs = get_body_position_dofs(sys, bd2)
    bd2_vel_dofs = get_body_velocity_dofs(sys, bd2)
    
    # Текущие углы и угловые скорости
    θ1 = state[bd1_pos_dofs[3]]
    θ2 = state[bd2_pos_dofs[3]]
    ω1 = state[bd1_vel_dofs[3]]
    ω2 = state[bd2_vel_dofs[3]]
    
    # Относительное смещение и скорость
    Δθ = θ1 - θ2 - spring.rest_angle
    Δω = ω1 - ω2
    
    # Полный момент (пружина + демпфер)
    τ = spring.stiffness * Δθ + spring.damping * Δω #+ spring.stiffness * Δθ^3
    
    # Добавляем в угловые ускорения (делим на инерцию)
    rhs[bd1_vel_dofs[3]] += -τ 
    rhs[bd2_vel_dofs[3]] += τ
end
# Вспомогательная функция для расчета момента пружины
function get_spring_moment(spring::TorsionalSpring, θ1::Float64, θ2::Float64)
    Δθ1 = θ1 - spring.rest_angle
    Δθ2 = θ2 - spring.rest_angle
    return spring.stiffness * (Δθ1 - Δθ2)
end

# Функция для получения энергии пружины
function get_spring_energy(spring::TorsionalSpring, θ1::Float64, θ2::Float64)
    Δθ1 = θ1 - spring.rest_angle
    Δθ2 = θ2 - spring.rest_angle
    Δθ_rel = Δθ1 - Δθ2
    return 0.5 * spring.stiffness * Δθ_rel^2
end
