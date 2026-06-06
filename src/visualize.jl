function Makie.lift(system, solution, body::Body2D, i::Observable)
    return lift(i) do value
        points = Vector{Point2f}(undef, 2)
        points[1:2] .= get_boundary_points(system, body, view(solution, :, value))
        return points
    end
end
function Makie.lift(system, solution, joint::AbstractJoint2D, i::Observable)
end


function Makie.lift(system, solution, joint::HingeJoint, i::Observable)
    p =  lift(i) do value
        point = get_hinge_point(system, joint, view(solution, :, value)) 
        return point;
    end
    return p;
end
function get_torsionalSpring_point(system::MBSystem2D, spring::Union{TorsionalSpring,LinearSpring}, state::AbstractVector{Float64})
    bd1 = spring.hinge.body1
    pos_dofs1 = get_body_position_dofs(system, bd1)
    _xi1 = state[pos_dofs1[1]]
    _yi1 = state[pos_dofs1[2]]
    _θi1 = state[pos_dofs1[3]]

    bd2 = spring.hinge.body2
    pos_dofs2 = get_body_position_dofs(system, bd2)
    _xi2 = state[pos_dofs2[1]]
    _yi2 = state[pos_dofs2[2]]
    _θi2 = state[pos_dofs2[3]]
    R = SA_F64[
        cos(_θi1) -sin(_θi1);
        sin(_θi1) cos(_θi1)
    ]
    v = R * spring.hinge.body1_hinge_point
    _xi = _xi1 + v[1]
    _yi = _yi1 + v[2]

    return Point2f(_xi ,_yi)
end

function get_Spring_point(system::MBSystem2D, spring::LinearSpring, state::AbstractVector{Float64})
    bd1 = spring.joint.body1
    pos_dofs1 = get_body_position_dofs(system, bd1)
    _xi1 = state[pos_dofs1[1]]
    _yi1 = state[pos_dofs1[2]]
    _θi1 = state[pos_dofs1[3]]

    bd2 = spring.joint.body2
    pos_dofs2 = get_body_position_dofs(system, bd2)
    _xi2 = state[pos_dofs2[1]]
    _yi2 = state[pos_dofs2[2]]
    _θi2 = state[pos_dofs2[3]]

    _xi = (_yi1 - _yi2 + tan(_θi2) * _xi2 - tan(_θi1) * _xi1) / (tan(_θi2) - tan(_θi1))
    _yi = tan(_θi2) * (_xi - _xi2) + _yi2

    return Point2f(_xi ,_yi)
end

function Makie.lift(system, solution, spring::LinearSpring, i::Observable)
    p = lift(i) do value
        point = get_Spring_point(system, spring, view(solution, :, value))
        bd1 = spring.joint.body1
        pos_dofs1 = get_body_position_dofs(system, bd1)
        _xi1, _yi1, _θi1 = view(solution, :, value)[pos_dofs1]

        bd2 = spring.joint.body2
        pos_dofs2 = get_body_position_dofs(system, bd2)
        _xi2, _yi2, _θi2 = view(solution, :, value)[pos_dofs2]

        t_spring = spring.vis_r / 2
        N = spring.vis_N
        x_range = LinRange(_xi1, _xi2, N)
        y_range = LinRange(_yi1, _yi2, N)
        points = Vector{Point2f}(undef, N)
        for j in 1:N
            if (j == 1)
                x = _xi1
                y = _yi1
            elseif (j == N)
                x = _xi2
                y = _yi2
            else
                x = x_range[j]
                y = y_range[j] + t_spring * (-1)^j
            end
            points[j] = Point2f(x, y)
        end
        return points;
    end
    return p
end

function Makie.lift(system, solution, spring::TorsionalSpring, i::Observable)
    p = lift(i) do value
        point = get_torsionalSpring_point(system, spring, view(solution, :, value)) 
        bd1 = spring.hinge.body1
        pos_dofs1 = get_body_position_dofs(system, bd1)
        _xi1, _yi1, _θi1 = view(solution, :, value)[pos_dofs1]

        bd2 = spring.hinge.body2
        pos_dofs2 = get_body_position_dofs(system, bd2)
        _xi2, _yi2, _θi2 = view(solution, :, value)[pos_dofs2]

        start_angel = _θi1 + π
        end_angel = _θi2 + 4*π

        r0 = spring.vis_r/2
        r1 = spring.vis_r
        N = 100

        t = LinRange(start_angel, end_angel, N)
        R = LinRange(r0, r1, N)
        x0 = point[1]
        y0 = point[2]

        points = Vector{Point2f}(undef, N)
        for j in 1:N
            x = R[j] * cos(t[j]) + x0  
            y = R[j] * sin(t[j]) + y0
            points[j] = Point2f(x,y)
        end
        return points;
    end
    return p
end

function draw!(ax, joint::AbstractConnector2D, system::MBSystem2D, solution, iter::Observable)
end
function draw!(ax, joint::HingeJoint, system::MBSystem2D, solution, iter::Observable)
    hinge_point = lift(system, solution, joint, iter);
    scatter!(ax, hinge_point);
end

function draw!(ax, joint::TorsionalSpring, system::MBSystem2D, solution, iter::Observable)
    hinge_point = lift(system, solution, joint, iter);
    lines!(ax, hinge_point);
end

function draw!(ax, joint::LinearSpring, system::MBSystem2D, solution, iter::Observable)
    hinge_point = lift(system, solution, joint, iter);
    lines!(ax, hinge_point);
end

# function draw!(ax, joint::FixedJoint, system::MBSystem2D, solution, iter::Observable)
#     fixed_point = lift(system, solution, joint, iter);
#     lines!(ax, fixed_point);
# end


function animate(sys::MBSystem2D, sol, time_span, filename; framerate=60, limits = (-1, 1, 1, 1))
    fig = Figure()
    iter = Observable(1)
    ax = Axis(fig[1, 1], aspect = DataAspect())

    for body in bodies(sys)
        bar = lift(sys, sol, body, iter)
        lines!(ax, bar)
    end

    for connector in connectors(sys)
        draw!(ax, connector, sys, sol, iter)
    end
    for force in sys.forces
        draw!(ax, force, sys, sol, iter)
    end

    limits!(ax, limits...)

    record(fig, filename, 1:5:length(time_span);
        framerate=framerate / 5) do t
        iter[] = t
    end
end


function draw_static(sys::MBSystem2D, sol;  limits = (-1, 1, 1, 1))
    fig = Figure()
    iter = Observable(1)
    ax = Axis(fig[1, 1], aspect = DataAspect())

    for body in bodies(sys)
        bar = lift(sys, sol, body, iter)
        lines!(ax, bar)
    end

    for connector in connectors(sys)
        draw!(ax, connector, sys, sol, iter)
    end
    for force in sys.forces
        draw!(ax, force, sys, sol, iter)
    end

    limits!(ax, limits...)

    return fig;
end