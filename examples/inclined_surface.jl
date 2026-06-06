using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using StaticArrays

# Parameters
const g = 9.81
angle_deg = 20.0
angle_rad = deg2rad(angle_deg)
direction = SA[cos(angle_rad), sin(angle_rad)]

# 1. Create bodies
# Fixed body (representing the track/ground)
fixed_body = Body2D(100.0, 10.0) 
# Sliding body
sliding_body = Body2D(1.0, 0.1)

# Apply gravity to sliding body
sliding_body.forces[2] = (state, t) -> -sliding_body.mass * g

# 2. Setup System
sys = MBSystem2D()
add!(sys, fixed_body)
add!(sys, sliding_body)

# 3. Create and configure joints
# Fix the first body in space
fixed_jnt = FixedJoint(fixed_body)
setposition!(fixed_jnt, SA[0.0, 0.0])
setrotation!(fixed_jnt, 0.0)

# Connect sliding body to fixed body via SliderJoint
slider_jnt = SliderJoint(fixed_body, sliding_body)
# Slider starts at origin of fixed body
set_position_on_first_body!(slider_jnt, SA[0.0, 0.0])
# Slider connected to center of sliding body
set_position_on_second_body!(slider_jnt, SA[0.0, 0.0])
# Set sliding axis direction (inclined at 20 degrees)
set_direction_on_first_body!(slider_jnt, direction)
set_direction_on_second_body!(slider_jnt, direction)

add!(sys, fixed_jnt)
add!(sys, slider_jnt)

# 4. Assemble and Simulate
if (!assemble!(sys))
    error("Assembling failed!")
end

# Initial state
# We need to satisfy constraints. 
# Body 1 is at (0,0,0). 
# Body 2 must be on the line passing through (0,0) with direction 'direction'.
initial_dist = 1.0
initial_pos = initial_dist .* direction
    
initial_state = zeros(number_of_dofs(sys))
# Body 2 position



# Simulation setup
func = sys.rhs
jacoby = (x) -> ForwardDiff.jacobian(func, x)

mass_matrix = get_mass_matrix(sys)

time_span = 0:0.02:2.0
sol = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol, initial_state, mass_matrix, func, jacoby, step(time_span))

# 5. Visualization
println("Simulated SliderJoint sliding under gravity at 20 degrees.")
animate(sys, sol, time_span, "out/slider_example.mp4"; framerate = 30, limits = (-10, 2, -10, 2))
