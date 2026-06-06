mutable struct Body2D <: AbstractBody2D
    mass::Float64
    inertia::Float64

    forces::Vector{Function}
    index::Int64
    length::Float64

    function Body2D(mass, inertia; length=1.0)
        return new(mass, inertia, [
                (x, t) -> 0,
                (x, t) -> 0,
                (x, t) -> 0
            ], -1, length)
    end
end
get_index(body::Body2D) = body.index
number_of_dofs(::Body2D) = 3
# get_lenght(body::Body2D)

function add!(sys::MBSystem2D, body::Body2D)
    push!(sys.bodies, body)
    # sys.bodiesnum += 1
    bodydofs = 2 * number_of_dofs(body)
    lbd = last_body_dof(sys)
    push!(sys.bodiesdofs, lbd + bodydofs)
    body.index = length(sys.bodies)
end

function get_body_position_dofs(sys::MBSystem2D, body::Body2D)
    bid = get_index(body)
    return SA[
        sys.bodiesdofs[bid]-5,
        sys.bodiesdofs[bid]-4,
        sys.bodiesdofs[bid]-3,
    ]

end
function get_body_velocity_dofs(sys::MBSystem2D, body::Body2D)
    bid = get_index(body)
    return SA[
        sys.bodiesdofs[bid]-2,
        sys.bodiesdofs[bid]-1,
        sys.bodiesdofs[bid]-0,
    ]
end

function get_body_generalized_dofs(sys::MBSystem2D, body::Body2D)
    bid = get_index(body)
    lg_dof = Int(sys.bodiesdofs[bid] ÷ 2)
    return SA[
        lg_dof - 2,
        lg_dof - 1,
        lg_dof - 0,
    ]
end

function add_to_rhs!(rhs, state, sys, body)
    last_dof = sys.bodiesdofs[body.index]

    position_dofs = [
        last_dof - 5,
        last_dof - 4,
        last_dof - 3,
    ]

    velocity_dofs = [
        last_dof - 2,
        last_dof - 1,
        last_dof,
    ]
    # assemble velocities
    rhs[position_dofs[1]] = state[velocity_dofs[1]]
    rhs[position_dofs[2]] = state[velocity_dofs[2]]
    rhs[position_dofs[3]] = state[velocity_dofs[3]]

    rhs[velocity_dofs[1]] += body.forces[1](state[[position_dofs; velocity_dofs]], state[end])
    rhs[velocity_dofs[2]] += body.forces[2](state[[position_dofs; velocity_dofs]], state[end])
    rhs[velocity_dofs[3]] += body.forces[3](state[[position_dofs; velocity_dofs]], state[end])
end

function get_boundary_points(sys::MBSystem2D, body::Body2D, dofs::AbstractVector{Float64})::Tuple{Point2f, Point2f}
    pos_inds = get_body_position_dofs(sys, body)
    x, y, θ = dofs[pos_inds]

    x1 = x - body.length/2 * cos(θ)
    x2 = x + body.length/2 * cos(θ)
    y1 = y - body.length/2 * sin(θ)
    y2 = y + body.length/2 * sin(θ)

    return Point2f(x1, y1), Point2f(x2, y2)
end
