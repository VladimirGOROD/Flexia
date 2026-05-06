
mutable struct FixedJoint <: AbstractJoint2D
    body::Body2D
    pos::SVector{2,Float64}
    θ::Float64
    index::Int64

    function FixedJoint(body::Body2D)
        return new(body, [0, 0], 0, -1)
    end
end

function setposition!(joint::FixedJoint, pos)
    joint.pos = SA[pos[1], pos[2]]
    return nothing
end

function setrotation!(joint::FixedJoint, θ)
    joint.θ = θ
    return nothing
end

number_of_dofs(::FixedJoint) = 3

function add!(sys::MBSystem2D, joint::AbstractJoint2D)
    push!(sys.joints, joint)
    # sys.jointsnum += 1;
    last_joint_dof = last_lm_dof(sys) + number_of_dofs(joint)
    push!(sys.lmdofs, last_joint_dof)
    setid!(joint, length(sys.joints))
end

function get_lms(sys::MBSystem2D, joint::FixedJoint)
    last_lm = sys.lmdofs[joint.index]
    return SA[
        last_lm-2,
        last_lm-1,
        last_lm-0,
    ]
end

function get_fixed_point(sys::MBSystem2D, joint::FixedJoint, state::AbstractVector{Float64})

    bd = joint.body
    pos_dofs = get_body_position_dofs(sys, bd)

    _xi = state[pos_dofs[1]]
    _yi = state[pos_dofs[2]]

    return Point2f(_xi, _yi)
end

function add_joint_to_rhs!(rhs, state, sys::MBSystem2D, joint::FixedJoint)
    body = joint.body
    last_body_dof = sys.bodiesdofs[body.index]

    position_dofs = [
        last_body_dof - 5,
        last_body_dof - 4,
        last_body_dof - 3,
    ]

    velocity_dofs = [
        last_body_dof - 2,
        last_body_dof - 1,
        last_body_dof,
    ]

    joint_dofs = get_lms(sys, joint)

    rhs[velocity_dofs[1]] += state[joint_dofs[1]]
    rhs[velocity_dofs[2]] += state[joint_dofs[2]]
    rhs[velocity_dofs[3]] += state[joint_dofs[3]]

    rhs[joint_dofs[1]] = state[position_dofs[1]] - joint.pos[1]
    rhs[joint_dofs[2]] = state[position_dofs[2]] - joint.pos[2]
    rhs[joint_dofs[3]] = state[position_dofs[3]] - joint.θ
end


mutable struct HingeJoint <: AbstractJoint2D
    body1::Body2D
    body1_hinge_point::SVector{2,Float64}
    body2::Body2D
    body2_hinge_point::SVector{2,Float64}
    index::Int64

    function HingeJoint(bd1::Body2D, bd2::Body2D)
        return new(bd1, SA[0.0, 0.0], bd2, SA[0.0, 0.0], -1)
    end
end

number_of_dofs(::HingeJoint) = 2
function setid!(jnt::T, index) where T<:AbstractJoint2D
    jnt.index = index
end
function set_position_on_first_body!(joint::HingeJoint, pos::SVector{2,Float64})
    joint.body1_hinge_point = pos
    return nothing
end

function set_position_on_second_body!(joint::HingeJoint, pos::SVector{2,Float64})
    joint.body2_hinge_point = pos
    return nothing
end

function get_lms(sys::MBSystem2D, joint::HingeJoint)
    last_lm = sys.lmdofs[joint.index]
    return SA[
        last_lm-1,
        last_lm-0,
    ]
end

function add_joint_to_rhs!(rhs, state, sys::MBSystem2D, joint::HingeJoint)
    bd1 = joint.body1
    bd2 = joint.body2

    bd1_p_dofs = get_body_position_dofs(sys, bd1)
    bd1_v_dofs = get_body_velocity_dofs(sys, bd1)
    bd2_p_dofs = get_body_position_dofs(sys, bd2)
    bd2_v_dofs = get_body_velocity_dofs(sys, bd2)

    _xi = state[bd1_p_dofs[1]]
    _yi = state[bd1_p_dofs[2]]
    _θi = state[bd1_p_dofs[3]]

    # vxi = state[bd1_v_dofs[1]]
    # vyi = state[bd1_v_dofs[2]]
    # _ωi = state[bd1_v_dofs[3]]

    _xj = state[bd2_p_dofs[1]]
    _yj = state[bd2_p_dofs[2]]
    _θj = state[bd2_p_dofs[3]]

    # vxj = state[bd2_v_dofs[1]]
    # vyj = state[bd2_v_dofs[2]]
    # _ωj = state[bd2_v_dofs[3]]

    lms = get_lms(sys, joint)
    λ1 = state[lms[1]]
    λ2 = state[lms[2]]

    xci = joint.body1_hinge_point[1]
    yci = joint.body1_hinge_point[2]
    xcj = joint.body2_hinge_point[1]
    ycj = joint.body2_hinge_point[2]

    rhs[bd1_v_dofs[1]] += λ1
    rhs[bd1_v_dofs[2]] += λ2
    rhs[bd1_v_dofs[3]] += (-xci * sin(_θi) - yci * cos(_θi)) * λ1 +
                          (xci * cos(_θi) - yci * sin(_θi)) * λ2

    rhs[bd2_v_dofs[1]] += -λ1
    rhs[bd2_v_dofs[2]] += -λ2
    rhs[bd2_v_dofs[3]] += (xcj * sin(_θj) + ycj * cos(_θj)) * λ1 +
                          (-xcj * cos(_θj) + ycj * sin(_θj)) * λ2

    rhs[lms[1]] = (_xi + xci * cos(_θi) - yci * sin(_θi)) -
                  (_xj + xcj * cos(_θj) - ycj * sin(_θj))
    rhs[lms[2]] = (_yi + xci * sin(_θi) + yci * cos(_θi)) -
                  (_yj + xcj * sin(_θj) + ycj * cos(_θj))


end

function get_hinge_point(system::MBSystem2D, joint::HingeJoint, state::AbstractVector{Float64})
    bd = joint.body1
    pos_dofs = get_body_position_dofs(system, bd)
    _xi = state[pos_dofs[1]]
    _yi = state[pos_dofs[2]]
    _θi = state[pos_dofs[3]]

    xci = joint.body1_hinge_point[1]
    yci = joint.body1_hinge_point[2]

    return Point2f(
        _xi + cos(_θi) * xci - sin(_θi) * yci,
        _yi + sin(_θi) * xci + cos(_θi) * yci
    )
end

mutable struct SliderJoint <: AbstractJoint2D
    body1::Body2D
    body1_position::SVector{2,Float64}
    body1_direction::SVector{2,Float64}
    body2::Body2D
    body2_position::SVector{2,Float64}
    body2_direction::SVector{2,Float64}
    alpha1::Float64
    alpha2::Float64
    index::Int64


    function SliderJoint(bd1::Body2D, bd2::Body2D)
        α1 = 0;
        α2 = 0;
        return new(bd1, SA[0.0, 0.0], SA[1.0, 0.0], bd2, SA[0.0, 0.0], SA[1.0, 0.0], α1, α2, -1)
    end
end

number_of_dofs(::SliderJoint) = 2

function get_lms(sys::MBSystem2D, joint::SliderJoint)
    last_lm = sys.lmdofs[joint.index]
    return SA[
        last_lm-1,
        last_lm-0,
    ]
end

function set_position_on_first_body!(joint::SliderJoint, pos::SVector{2,Float64})
    joint.body1_position = pos
    return nothing
end

function set_position_on_second_body!(joint::SliderJoint, pos::SVector{2,Float64})
    joint.body2_position = pos
    return nothing
end

function set_direction_on_first_body!(joint::SliderJoint, dir::SVector{2,Float64})
    joint.alpha1 = atan(dir[2], dir[1])
    joint.body1_direction = dir
    return nothing
end

function set_direction_on_second_body!(joint::SliderJoint, dir::SVector{2,Float64})
    joint.alpha2 = atan(dir[2], dir[1])
    joint.body2_direction = dir
    return nothing
end

function add_joint_to_rhs!(rhs, state, sys::MBSystem2D, joint::SliderJoint)
    bd1 = joint.body1
    bd2 = joint.body2

    bd1_p_dofs = get_body_position_dofs(sys, bd1)
    bd1_v_dofs = get_body_velocity_dofs(sys, bd1)
    bd2_p_dofs = get_body_position_dofs(sys, bd2)
    bd2_v_dofs = get_body_velocity_dofs(sys, bd2)

    _xi = state[bd1_p_dofs[1]]
    _yi = state[bd1_p_dofs[2]]
    _θi = state[bd1_p_dofs[3]]

    # vxi = state[bd1_v_dofs[1]]
    # vyi = state[bd1_v_dofs[2]]
    # ωi = state[bd1_v_dofs[3]]

    _xj = state[bd2_p_dofs[1]]
    _yj = state[bd2_p_dofs[2]]
    _θj = state[bd2_p_dofs[3]]

    # vxj = state[bd2_v_dofs[1]]
    # vyj = state[bd2_v_dofs[2]]
    # ωj = state[bd2_v_dofs[3]]

    lms = get_lms(sys, joint)
    λ1 = state[lms[1]]
    λ2 = state[lms[2]]

    xci = joint.body1_position[1]
    yci = joint.body1_position[2]
    xdi = joint.body1_direction[1]
    ydi = joint.body1_direction[2]

    xcj = joint.body2_position[1]
    ycj = joint.body2_position[2]
    xdj = joint.body2_direction[1]
    ydj = joint.body2_direction[2]

    αi = joint.alpha1
    αj = joint.alpha2
    
    # Constraint equations for slider joint:
    # 1. Perpendicular distance between points is zero
    # 2. Rotation difference (directions aligned)

    # Position of hinge points in global coordinates
    xpi = _xi + xci * cos(_θi) - yci * sin(_θi)
    ypi = _yi + xci * sin(_θi) + yci * cos(_θi)

    xpj = _xj + xcj * cos(_θj) - ycj * sin(_θj)
    ypj = _yj + xcj * sin(_θj) + ycj * cos(_θj)


    # normal to direction in terms of first body in global CS
    G_xni = -xdi*sin(_θi) - ydi*cos(_θi)
    G_yni =  xdi*cos(_θi) - ydi*sin(_θi)

    # Constraint 1: perpendicular distance
    rhs[lms[1]] = (xpj - xpi) * G_xni + (ypj - ypi) * G_yni

    # Constraint 2: direction alignment (rotation difference)
    # dx_gj * dy_gi - dy_gj * dx_gi = 0 (cross product)
    rhs[lms[2]] = _θj + αj - (_θi + αi)

    # Velocity constraints
    # For constraint 1: d/dt of the perpendicular distance
    rhs[bd1_v_dofs[1]] += λ1 * G_xni 
    rhs[bd1_v_dofs[2]] += λ1 * G_yni
    rhs[bd1_v_dofs[3]] += λ1 * ((_yi - ypj) * G_xni + (xpj - _xi) * G_yni)

    rhs[bd2_v_dofs[1]] += -λ1 * G_xni
    rhs[bd2_v_dofs[2]] += -λ1 * G_yni
    rhs[bd2_v_dofs[3]] += -λ1 * ((ypj - _yj ) * G_xni + (_xj - xpj) * G_yni)

    # For constraint 2: d/dt of direction alignment
    rhs[bd1_v_dofs[3]] += -λ2 
    rhs[bd2_v_dofs[3]] += λ2
end