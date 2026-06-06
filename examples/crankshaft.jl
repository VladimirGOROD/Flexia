using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using StaticArrays


# Parameters
const g = 9.81

# Crankshaft mechanism example demonstrating the use of HingeJoint and SliderJoint
#
# This example models a crankshaft mechanism with 4 bodies:
# 1. Ground: Fixed at (0,0), length=0.01
# 2. Crank: Rotating at (0.5,0), length=1.0
# 3. Horizontal slider: Sliding at (3,0), length=0.25
# 4. Connection rod: Connecting crank to slider at (2,0), length=2.0
#
# Joints:
# 1. Ground-to-crank: HingeJoint at point (0,0) on ground, (-0.5,0) on crank
# 2. Crank-to-connector: HingeJoint at point (0.5,0) on crank, (-1,0) on connector
# 3. Connector-to-slider: HingeJoint at point (1,0) on connector, (0,0) on slider
# 4. Slider-to-ground: SliderJoint at point (0,0) on slider, (3,0) on ground, axis (1,0)

# Create bodies
ground = Body2D(1e6, 1e6; length=0.01)  # Large mass/inertia to simulate fixed ground
crank = Body2D(10.0, 1.0; length=1.0)
slider = Body2D(5.0, 0.1; length=0.25)
connector = Body2D(5.0, 0.5; length=2.0)

# Apply gravity
crank.forces[2] = (x, t) -> -crank.mass * g
slider.forces[2] = (x, t) -> -slider.mass * g
slider.forces[1] = (x, t) -> 2050*sin(4*t)
connector.forces[2] = (x, t) -> -connector.mass * g

# Create system
sys = MBSystem2D()
add!(sys, ground)
add!(sys, crank)
add!(sys, connector)
add!(sys, slider)

# Create joints
# 1. Fix ground
ground_joint = FixedJoint(ground)
setposition!(ground_joint, SA[0.0, 0.0])
setrotation!(ground_joint, 0.0)

# 2. Ground to crank hinge
crank_ground_hinge = HingeJoint(ground, crank)
set_position_on_first_body!(crank_ground_hinge, SA[0.0, 0.0])
set_position_on_second_body!(crank_ground_hinge, SA[-0.5, 0.0])

# 3. Crank to connector hinge
crank_connector_hinge = HingeJoint(crank, connector)
set_position_on_first_body!(crank_connector_hinge, SA[0.5, 0.0])
set_position_on_second_body!(crank_connector_hinge, SA[-1.0, 0.0])

# 4. Connector to slider hinge
connector_slider_hinge = HingeJoint(connector, slider)
set_position_on_first_body!(connector_slider_hinge, SA[1.0, 0.0])
set_position_on_second_body!(connector_slider_hinge, SA[0.0, 0.0])

# 5. Slider to ground slider joint
slider_ground_slider = SliderJoint(ground, slider)
set_position_on_first_body!(slider_ground_slider, SA[3.0, 0.0])
set_position_on_second_body!(slider_ground_slider, SA[0.0, 0.0])
set_direction_on_first_body!(slider_ground_slider, SA[1.0, 0.0])
set_direction_on_second_body!(slider_ground_slider, SA[1.0, 0.0])

# Add joints to system
add!(sys, ground_joint)
add!(sys, crank_ground_hinge)
add!(sys, crank_connector_hinge)
add!(sys, connector_slider_hinge)
add!(sys, slider_ground_slider)

# Assemble system
assemble!(sys)
println("System assembled successfully")

# Set up simulation
func = sys.rhs
jacoby = (x) -> ForwardDiff.jacobian(func, x)

# Initial conditions - need to satisfy constraints
initial_state = zeros(number_of_dofs(sys))

# Set initial positions
# Ground: already fixed at (0,0,0)
# Crank: center at (0.5, 0), rotation 0
crank_x, crank_y, crank_θ = get_body_position_dofs(sys, crank)
initial_state[crank_x] = 0.5
initial_state[crank_y] = 0.0
initial_state[crank_θ] = 0.0

# Connector: center at (2,0), rotation such that it connects properly
connector_x, connector_y, connector_θ = get_body_position_dofs(sys, connector)
initial_state[connector_x] = 2.0
initial_state[connector_y] = 0.0
initial_state[connector_θ] = 0.0

# Slider: center at (3,0), rotation 0 (aligned with x-axis)
slider_x, slider_y, slider_θ = get_body_position_dofs(sys, slider)
initial_state[slider_x] = 3.0
initial_state[slider_y] = 0.0
initial_state[slider_θ] = 0.0

# Mass matrix
# mass_matrix = zeros(number_of_dofs(sys), number_of_dofs(sys))
# for i in 1:last_body_dof(sys)
    # mass_matrix[i, i] = 1.0
# end
mass_matrix = sys.mass
# Simulation parameters
time_span = 0:0.01:10.0
sol = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))

# Solve
cros!(sol, initial_state, mass_matrix, func, jacoby, step(time_span))

# Animate
animate(sys, sol, time_span, "out/crankshaft.mp4"; framerate=30, limits=(-2, 5, -2, 2))
