mutable struct PositionMotor2D <: AbstractPositionActuator2D
    joint::HingeJoint
    target_angle::Float64
    index::Int64

    function PositionMotor2D(joint::HingeJoint, target_angle::Float64)
        return new(joint, target_angle, -1)
    end
end



function settarget!(act::PositionMotor2D, angle::Float64)
    act.target_angle = angle
    return nothing
end

number_of_dofs(::PositionMotor2D) = 1

function get_lms(sys::MBSystem2D, act::PositionMotor2D)
    
    last_lm = sys.lmdofs[act.index]
    return SA[last_lm]
end

function add_to_rhs!(rhs, state, sys::MBSystem2D, act::PositionMotor2D)
    joint = act.joint
    bd1 = joint.body1
    bd2 = joint.body2

    bd1_p_dofs = get_body_position_dofs(sys, bd1)
    bd1_v_dofs = get_body_velocity_dofs(sys, bd1)
    bd2_p_dofs = get_body_position_dofs(sys, bd2)
    bd2_v_dofs = get_body_velocity_dofs(sys, bd2)

    θ1 = state[bd1_p_dofs[3]]
    θ2 = state[bd2_p_dofs[3]]

    ω1 = state[bd1_v_dofs[3]]
    ω2 = state[bd2_v_dofs[3]]

    lms = get_lms(sys, act)
    λ = state[lms[1]]


    # Constraint: θ2 - θ1 - target = 0
    rhs[lms[1]] = θ2 - θ1

    rhs[bd1_p_dofs[3]] -= λ
    rhs[bd2_p_dofs[3]] += λ
end

function propagate_targets!(sys::MBSystem2D, act::PositionMotor2D)
    lms = get_lms(sys, act)
    sys.targets[lms[1]] = act.target_angle
    return nothing
end


mutable struct PositionLinearActuator2D <: AbstractPositionActuator2D
    joint::SliderJoint
    target_position::Float64
    index::Int64

    function PositionLinearActuator2D(joint::SliderJoint, target_angle::Float64)
        return new(joint, target_angle, -1)
    end
end



function settarget!(act::PositionLinearActuator2D, position::Float64)
    act.target_position = position
    return nothing
end

number_of_dofs(::PositionLinearActuator2D) = 1

function get_lms(sys::MBSystem2D, act::PositionLinearActuator2D)
    
    last_lm = sys.lmdofs[act.index]
    return SA[last_lm]
end


function propagate_targets!(sys::MBSystem2D, act::PositionLinearActuator2D)
    lms = get_lms(sys, act)
    sys.targets[lms[1]] = act.target_position
    return nothing
end

function add_to_rhs!(rhs, state, sys::MBSystem2D, act::PositionLinearActuator2D)
    joint = act.joint
    bd1 = joint.body1
    bd2 = joint.body2

    bd1_p_dofs = get_body_position_dofs(sys, bd1)
    bd2_p_dofs = get_body_position_dofs(sys, bd2)

    _xi = state[bd1_p_dofs[1]]
    _yi = state[bd1_p_dofs[2]]
    _θi = state[bd1_p_dofs[3]]

    _xj = state[bd2_p_dofs[1]]
    _yj = state[bd2_p_dofs[2]]
    _θj = state[bd2_p_dofs[3]]

    xci = joint.body1_position[1]
    yci = joint.body1_position[2]
    xdi = joint.body1_direction[1]
    ydi = joint.body1_direction[2]

    xcj = joint.body2_position[1]
    ycj = joint.body2_position[2]

    αi = joint.alpha1
    αj = joint.alpha2
    
    xpi = _xi + xci * cos(_θi) - yci * sin(_θi)
    ypi = _yi + xci * sin(_θi) + yci * cos(_θi)

    xpj = _xj + xcj * cos(_θj) - ycj * sin(_θj)
    ypj = _yj + xcj * sin(_θj) + ycj * cos(_θj)

    # direction in terms of first body in global CS
    N_xni = xdi*cos(_θi) - ydi*sin(_θi)
    N_yni = xdi*sin(_θi) + ydi*cos(_θi)

    lms = get_lms(sys, act)
    λ = state[lms[1]]

    rhs[lms[1]] = (xpj-xpi) * N_xni + (ypj-ypi) * N_yni

    rhs[bd1_p_dofs[1]] += λ * (-N_xni)
    rhs[bd1_p_dofs[2]] += λ * (-N_yni)
    rhs[bd1_p_dofs[3]] += λ * ((ypj - _yi)*(N_xni) - (xpj - _xi) * N_yni)
    
    rhs[bd2_p_dofs[1]] += λ * (N_xni)
    rhs[bd2_p_dofs[2]] += λ * (N_yni)
    rhs[bd2_p_dofs[3]] += λ * ((_yj - ypj)*N_xni + (_xj - xpj)*N_yni)
end