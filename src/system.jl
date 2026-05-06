mutable struct MBSystem2D
    bodies::Vector{AbstractBody2D}
    joints::Vector{AbstractJoint2D}
    forces::Vector{AbstractForce2D}
    # bodiesnum::Int64
    # jointsnum::Int64
    bodiesdofs::Vector{Int64}
    lmdofs::Vector{Int64}
    assembled::Bool
    rhs::Function
    jacobian::Function

    function MBSystem2D()
        default_rhs = (x) -> nothing
        default_jacobian = (x) -> nothing
        return new([], [], [], [], [], false, default_rhs, default_jacobian)
    end
end
number_of_bodies(sys::MBSystem2D) = length(sys.bodies)
bodies(sys::MBSystem2D) = sys.bodies
joints(sys::MBSystem2D) = sys.joints
function last_body_dof(sys::MBSystem2D)
    if (isempty(sys.bodiesdofs))
        return 0
    end
    return (sys.bodiesdofs)[end]
end
function last_lm_dof(sys::MBSystem2D)
    if (isempty(sys.lmdofs))
        return last_body_dof(sys)
    else
        return sys.lmdofs[end]
    end
end

number_of_dofs(sys) = last_lm_dof(sys)

function assemble!(sys)
    state_length = number_of_dofs(sys)


    sys.rhs = (state) -> begin
        ret = similar(state)
        fill!(ret, zero(state[1]))
        for body in sys.bodies
            add_body_to_rhs!(ret, state, sys, body)
        end
        for joint in sys.joints
            add_joint_to_rhs!(ret, state, sys, joint)
        end
        for force in sys.forces
            add_force_to_rhs!(ret, state, sys, force)
        end
        return ret
    end

    sys.jacobian = (state) -> begin
        ret = zeros(state_length, state_length)
        for (i, body) in enumerate(sys.bodies)
            last_dof = sys.bodiesdofs[i]

            position_dofs = [
                last_dof - 5,
                last_dof - 4,
                last_dof - 3,
            ]
            for p in position_dofs
                ret[p, p] = 1
            end
            velocity_dofs = [
                last_dof - 2,
                last_dof - 1,
                last_dof,
            ]
            for (j, v) in enumerate(velocity_dofs)
                ret[v, last_dof-5:last_dof] = ForwardDiff().gradient(
                    body.forces[j],
                    state[last_dof-5:last_dof])
                ret[v, velocity_dofs] = ForwardDiff().gradient(
                    (x) -> body.forces[j](state[position_dofs], x),
                    state[velocity_dofs])
            end
        end
        return ret
    end

    sys.assembled = true
end


# function draw!(ax, sys::MBSystem2D, state::Vector)
#     for body in bodies(sys)
#         rod = get_boundary_points(sys, bd1, state)
#         draw!(ax, sys, body, state);
#     end
#     for joint in joints(sys)
#         draw!()
#     end
# end