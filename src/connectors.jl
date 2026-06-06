function setid!(conn::AbstractConnector2D, index)
    conn.index = index
end

function add!(sys::MBSystem2D, connector::CT) where CT<:AbstractConnector2D
    push!(sys.connectors, connector)
    # sys.connectorsnum += 1;
    last_connector_dof = last_lm_dof(sys) + number_of_dofs(connector)
    push!(sys.lmdofs, last_connector_dof)
    setid!(connector, length(sys.connectors))
end

function propagate_targets!(sys::MBSystem2D, connector::T) where T<:AbstractConnector2D
    return nothing
end

function compute_kinematic_residual!(
    residual::Vector{Float64}, 
    generalized_coordinates::Vector{Float64}, 
    sys::MBSystem2D, 
    connector::AbstractConnector2D)
    return nothing
end