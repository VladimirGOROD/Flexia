mutable struct TorsionalSpring <: AbstractForce2D
    hinge::HingeJoint
    stiffness::Float64
    rest_angle::Float64
    damping::Float64
    vis_r::Float64
    
    function TorsionalSpring(hinge::HingeJoint, stiffness::Float64=0.0, rest_angle::Float64=0.0, damping::Float64=0.0, vis_r =0.3)
        return new(hinge, stiffness, rest_angle, damping, vis_r)
    end

    # function TorsionalSpring(body1::Body2D, body2::Body2D, stiffness::Float64=1.0, rest_angle::Float64=0.0, damping::Float64=0.0, vis_r =0.3)
    #     return new(body1, body2, stiffness, rest_angle, damping, vis_r)
    # end
end
function TorsionalSpring(hinge::HingeJoint; stiffness::Float64=0.0, rest_angle::Float64=0.0, damping::Float64=0.0, vis_r =0.3)
    return TorsionalSpring(hinge, stiffness, rest_angle, damping, vis_r)
end
number_of_dofs(::TorsionalSpring) = 0  # пружина не добавляет лагранжевых множителей

function add!(sys::MBSystem2D, spring::TorsionalSpring)
    push!(sys.forces, spring)
end

function add_to_rhs!(rhs, state, sys::MBSystem2D, spring::TorsionalSpring)
    

    bd1 = spring.hinge.body1
    bd2 = spring.hinge.body2
    
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
    
    # Добавляем в угловые ускорения
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

mutable struct LinearSpring <: AbstractForce2D
    joint::SliderJoint
    stiffness::Float64
    x_dist::Float64
    damping::Float64
    vis_r::Float64
    vis_N::Int
    
    function LinearSpring(joint::SliderJoint, stiffness::Float64=1.0, x_dist::Float64=0.0, damping::Float64=0.0, vis_r =0.3, vis_N = 4)
        return new(joint, stiffness, x_dist, damping, vis_r, vis_N)
    end
end

function LinearSpring(joint::SliderJoint; stiffness::Float64=1.0, x_dist::Float64=0.0, damping::Float64=0.0, vis_r =0.3, vis_N = 4)
    return LinearSpring(joint, stiffness, x_dist, damping, vis_r, vis_N)
end

number_of_dofs(::LinearSpring) = 0  # пружина не добавляет лагранжевых множителей

function add!(sys::MBSystem2D, spring::LinearSpring)
    push!(sys.forces, spring)
end

function add_to_rhs!(rhs, state, sys::MBSystem2D, spring::LinearSpring)


    joint = spring.joint
    bd1 = joint.body1
    bd2 = joint.body2

    bd1_p_dofs = get_body_position_dofs(sys, bd1)
    bd2_p_dofs = get_body_position_dofs(sys, bd2)

    bd1_v_dofs = get_body_velocity_dofs(sys, bd1)
    bd2_v_dofs = get_body_velocity_dofs(sys, bd2)

    _xi = state[bd1_p_dofs[1]]
    _yi = state[bd1_p_dofs[2]]
    _θi = state[bd1_p_dofs[3]]

    _xj = state[bd2_p_dofs[1]]
    _yj = state[bd2_p_dofs[2]]
    _θj = state[bd2_p_dofs[3]]


    vxi = state[bd1_p_dofs[1]]
    vyi = state[bd1_p_dofs[2]]

    vxj = state[bd2_p_dofs[1]]
    vyj = state[bd2_p_dofs[2]]


    xci = joint.body1_position[1]
    yci = joint.body1_position[2]
    xdi = joint.body1_direction[1]
    ydi = joint.body1_direction[2]

    xcj = joint.body2_position[1]
    ycj = joint.body2_position[2]

    
    xpi = _xi + xci * cos(_θi) - yci * sin(_θi)
    ypi = _yi + xci * sin(_θi) + yci * cos(_θi)

    xpj = _xj + xcj * cos(_θj) - ycj * sin(_θj)
    ypj = _yj + xcj * sin(_θj) + ycj * cos(_θj)

    # direction in terms of first body in global CS
    N_xni = xdi*cos(_θi) - ydi*sin(_θi)
    N_yni = xdi*sin(_θi) + ydi*cos(_θi)

    Δx = (xpj - xpi) * N_xni + (ypj - ypi) * N_yni - spring.x_dist
    Δv = (vxj - vxi) * N_xni + (vyj - vyi) * N_yni

    F = spring.stiffness * Δx #- spring.damping * Δv #+ spring.stiffness * Δθ^3
    Fx = F * N_xni
    Fy = F * N_yni

    rhs[bd1_v_dofs[1]] += Fx
    rhs[bd2_v_dofs[1]] += -Fx

    rhs[bd1_v_dofs[2]] += Fy
    rhs[bd2_v_dofs[2]] += -Fy
end