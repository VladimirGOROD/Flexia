using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using StaticArrays

ground = Body2D(1e6, 1e6; length=0.01)  # Large mass/inertia to simulate fixed ground
const g = 9.81
crank_len = 1.0;
conn_len = 1.0;
table_heigh = 1.5;

H = 1.5;
L = 3.0
crank1 = Body2D(1, 1; length=crank_len)
crank2 = Body2D(1, 1; length=crank_len)
crank3 = Body2D(1, 1; length=crank_len)
crank4 = Body2D(1, 1.001; length=crank_len)

connector1 = Body2D(1, 1; length=conn_len)
connector2 = Body2D(1, 1; length=conn_len)
connector3 = Body2D(1, 1; length=conn_len)
connector4 = Body2D(1, 1; length=conn_len)

table1 = Body2D(1, 1; length=table_heigh)
table2 = Body2D(1, 1; length=table_heigh)

crank1.forces[2] = (x, t) -> -crank1.mass * g
crank2.forces[2] = (x, t) -> -crank2.mass * g
crank3.forces[2] = (x, t) -> -crank3.mass * g
crank4.forces[2] = (x, t) -> -crank3.mass * g

connector1.forces[2] = (x, t) -> -connector1.mass * g
connector2.forces[2] = (x, t) -> -connector2.mass * g
connector3.forces[2] = (x, t) -> -connector4.mass * g
connector4.forces[2] = (x, t) -> -connector4.mass * g
table1.forces[2] = (x, t) -> -table1.mass * g
table2.forces[2] = (x, t) -> -table2.mass * g



ground_joint = FixedJoint(ground)
setposition!(ground_joint, SA[0.0, 0.0])
setrotation!(ground_joint, π)

hinge_gr2cr1 = HingeJoint(ground, crank1)
set_position_on_first_body!(hinge_gr2cr1, SA[-L/2, H/2])
set_position_on_second_body!(hinge_gr2cr1, SA[-crank_len/2, 0])

hinge_gr2cr2 = HingeJoint(ground, crank2)
set_position_on_first_body!(hinge_gr2cr2, SA[-L/2, -H/2])
set_position_on_second_body!(hinge_gr2cr2, SA[-crank_len/2, 0])

hinge_gr2cr3 = HingeJoint(ground, crank3)
set_position_on_first_body!(hinge_gr2cr3, SA[L/2, -H/2])
set_position_on_second_body!(hinge_gr2cr3, SA[crank_len/2, 0])

hinge_gr2cr4 = HingeJoint(ground, crank4)
set_position_on_first_body!(hinge_gr2cr4, SA[L/2, H/2])
set_position_on_second_body!(hinge_gr2cr4, SA[crank_len/2, 0])

hinge_cr2con1 = HingeJoint(crank1, connector1)
set_position_on_first_body!(hinge_cr2con1, SA[crank_len/2, 0])
set_position_on_second_body!(hinge_cr2con1, SA[-conn_len/2, 0])

hinge_cr2con2 = HingeJoint(crank2, connector2)
set_position_on_first_body!(hinge_cr2con2, SA[crank_len/2, 0])
set_position_on_second_body!(hinge_cr2con2, SA[-conn_len/2, 0])

hinge_cr2con3 = HingeJoint(crank3, connector3)
set_position_on_first_body!(hinge_cr2con3, SA[-crank_len/2, 0])
set_position_on_second_body!(hinge_cr2con3, SA[conn_len/2, 0])

hinge_cr2con4 = HingeJoint(crank4, connector4)
set_position_on_first_body!(hinge_cr2con4, SA[-crank_len/2, 0])
set_position_on_second_body!(hinge_cr2con4, SA[conn_len/2, 0])


hinge_con2tab1 = HingeJoint(connector1, table1)
set_position_on_first_body!(hinge_con2tab1, SA[conn_len/2, 0])
set_position_on_second_body!(hinge_con2tab1, SA[table_heigh/2, 0])

hinge_con2tab2 = HingeJoint(connector2, table1)
set_position_on_first_body!(hinge_con2tab2, SA[conn_len/2, 0])
set_position_on_second_body!(hinge_con2tab2, SA[-table_heigh/2, 0])

hinge_con2tab3 = HingeJoint(connector3, table2)
set_position_on_first_body!(hinge_con2tab3, SA[-conn_len/2, 0])
set_position_on_second_body!(hinge_con2tab3, SA[-table_heigh/2, 0])

hinge_con2tab4 = HingeJoint(connector4, table2)
set_position_on_first_body!(hinge_con2tab4, SA[-conn_len/2, 0])
set_position_on_second_body!(hinge_con2tab4, SA[table_heigh/2, 0])

slider = SliderJoint(table1, table2)
set_position_on_first_body!(slider, SA[0.0, 0.0])
set_position_on_second_body!(slider, SA[0.0, 0.0])
set_direction_on_first_body!(slider, SA[0.0, -1.0])
set_direction_on_second_body!(slider, SA[0.0, -1.0])

sys = MBSystem2D()
add!(sys, ground)
add!(sys, crank1)
add!(sys, crank2)
add!(sys, crank3)
add!(sys, crank4)
add!(sys, connector1)
add!(sys, connector2)
add!(sys, connector3)
add!(sys, connector4)
add!(sys, table1)
add!(sys, table2)

add!(sys, hinge_gr2cr1)
add!(sys, hinge_gr2cr2)
add!(sys, hinge_gr2cr3)
add!(sys, hinge_gr2cr4)
add!(sys, hinge_cr2con1)
add!(sys, hinge_cr2con2)
add!(sys, hinge_cr2con3)
add!(sys, hinge_cr2con4)
add!(sys, hinge_con2tab1)
add!(sys, hinge_con2tab2)
add!(sys, hinge_con2tab3)
add!(sys, hinge_con2tab4)
add!(sys, slider)

assemble!(sys)
println("System assembled successfully")

func = sys.rhs
jacoby = (x) -> ForwardDiff.jacobian(func, x)
sys.jacobian = jacoby

# Compute initial positions
initial_state = zeros(number_of_dofs(sys))
func(initial_state)
# Ground
ground_x, ground_y, ground_θ = get_body_position_dofs(sys, ground)
initial_state[ground_x] = 0.0
initial_state[ground_y] = 0.0
initial_state[ground_θ] = 0.0

# Crank angles from motors
θ_crank1 = π/3
θ_crank2 = π/3
θ_crank3 = -π/3
θ_crank4 = -π/3
# Crank1
crank1_x, crank1_y, crank1_θ = get_body_position_dofs(sys, crank1)
initial_state[crank1_x] = -L/2 + (crank_len/2) * cos(θ_crank1)
initial_state[crank1_y] = H/2 + (crank_len/2) * sin(θ_crank1)
initial_state[crank1_θ] = θ_crank1

# Crank2
crank2_x, crank2_y, crank2_θ = get_body_position_dofs(sys, crank2)
initial_state[crank2_x] = -L/2 + (crank_len/2) * cos(θ_crank2)
initial_state[crank2_y] = -H/2 + (crank_len/2) * sin(θ_crank2)
initial_state[crank2_θ] = θ_crank2

# Crank3
crank3_x, crank3_y, crank3_θ = get_body_position_dofs(sys, crank3)
initial_state[crank3_x] = L/2 - (crank_len/2) * cos(θ_crank3)
initial_state[crank3_y] = -H/2 - (crank_len/2) * sin(θ_crank3)
initial_state[crank3_θ] = θ_crank3

# Crank4
crank4_x, crank4_y, crank4_θ = get_body_position_dofs(sys, crank4)
initial_state[crank4_x] = L/2 - (crank_len/2) * cos(θ_crank4)
initial_state[crank4_y] = H/2 - (crank_len/2) * sin(θ_crank4)
initial_state[crank4_θ] = θ_crank4

# Connectors
θ_conn1 = -θ_crank1
θ_conn2 = -θ_crank2
θ_conn3 = -θ_crank3
θ_conn4 = -θ_crank4

connector1_x, connector1_y, connector1_θ = get_body_position_dofs(sys, connector1)
initial_state[connector1_x] = initial_state[crank1_x] + cos(θ_crank1) * crank_len/2 + cos(θ_conn1) * (conn_len/2)
initial_state[connector1_y] = initial_state[crank1_y] + sin(θ_crank1) * crank_len/2 + sin(θ_conn1) * (conn_len/2)
initial_state[connector1_θ] = θ_conn1

connector2_x, connector2_y, connector2_θ = get_body_position_dofs(sys, connector2)
initial_state[connector2_x] = initial_state[crank2_x] + cos(θ_crank2) * crank_len/2 + cos(θ_conn2) * (conn_len/2)
initial_state[connector2_y] = initial_state[crank2_y] + sin(θ_crank2) * crank_len/2 + sin(θ_conn2) * (conn_len/2)
initial_state[connector2_θ] = θ_conn2

connector3_x, connector3_y, connector3_θ = get_body_position_dofs(sys, connector3)
initial_state[connector3_x] = initial_state[crank3_x] - cos(θ_crank3) * (crank_len/2) - cos(θ_conn3) * (conn_len/2)
initial_state[connector3_y] = initial_state[crank3_y] - sin(θ_crank3) * (crank_len/2) - sin(θ_conn3) * (conn_len/2)
initial_state[connector3_θ] = θ_conn3

connector4_x, connector4_y, connector4_θ = get_body_position_dofs(sys, connector4)
initial_state[connector4_x] = initial_state[crank4_x] - cos(θ_crank4) * (crank_len/2) - cos(θ_conn4) * (conn_len/2)
initial_state[connector4_y] = initial_state[crank4_y] - sin(θ_crank4) * (crank_len/2) - sin(θ_conn4) * (conn_len/2)
initial_state[connector4_θ] = θ_conn4

# Tables
θ_table1 = π/2
θ_table2 = π/2

table1_x, table1_y, table1_θ = get_body_position_dofs(sys, table1)
# Hinge from connector1
conn1_hinge_x = initial_state[connector1_x] + cos(θ_conn1) * (conn_len/2)
conn1_hinge_y = initial_state[connector1_y] + sin(θ_conn1) * (conn_len/2)
initial_state[table1_x] = conn1_hinge_x
initial_state[table1_y] = conn1_hinge_y - table_heigh/2
initial_state[table1_θ] = θ_table1

table2_x, table2_y, table2_θ = get_body_position_dofs(sys, table2)
conn3_hinge_x = initial_state[connector3_x] - cos(θ_conn3) * (conn_len/2)
conn3_hinge_y = initial_state[connector3_y] - sin(θ_conn3) * (conn_len/2)
initial_state[table2_x] = conn3_hinge_x
initial_state[table2_y] = conn3_hinge_y + table_heigh/2
initial_state[table2_θ] = θ_table2

println("Initial positions computed")

# Simulation parameters
time_span = 0:0.01:5

# Solve
sol = simulate(sys, initial_state, time_span)

# Animate
animate(sys, sol, time_span, "out/parmech.mp4"; framerate= floor(Int, 1.0 /step(time_span)), limits=(-4, 5, -2, 2))


