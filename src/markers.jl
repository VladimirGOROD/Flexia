mutable struct Marker <: AbstractMarker2D
    marker_name::String
    body_index::Int64
    point_on_body::AbstractVector{Float64}
    
    function first_marker(name::String, body::Body2D)
        marker_name = name
        body_index = get_index(body::Body2D)
        point_on_body = get_boundary_points(sys::MBSystem2D, body::Body2D, dofs::AbstractVector{Float64})
        return new(name, body_index)
    end

    function second_marker(name::String, body::Body2D)
        marker_name = name
        body_index = get_index(body)
        point_on_body = get_boundary_points(sys::MBSystem2D, body::Body2D, dofs::AbstractVector{Float64})
        return new(name, body_index)
    end
end