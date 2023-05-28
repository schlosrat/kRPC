import krpc
import time
import math
import datetime

## define functions we may need to call
# Define some vector (tuple) math functions to make things easier
def vector_subtract(u, v): # subtract vector v from vector u
	return tuple(map(lambda i, j: i - j, u, v))

def vector_add(u, v): # add vector v and vector u
	return tuple(map(lambda i, j: i + j, u, v))

def vector_mult_s(u, s): # multiply vector u by scalar s
	return tuple(x*s for x in u)

def vector_div_s(u, s): # divide vector u by scalar s
	return tuple(x/s for x in u)

def magnitude(u): # compute magnitude of vector u
	return math.sqrt(sum(i*i for i in u))

# def magnitude(v):
# 	return math.sqrt(dot_product(v, v))

def cross_product(u, v): # compute vector cross product of u x v
	return (u[1]*v[2] - u[2]*v[1],
	u[2]*v[0] - u[0]*v[2],
	u[0]*v[1] - u[1]*v[0])

# def dot_product(u, v): # Compute scalar dot product of u and v
# 	return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def dot_product(u, v):
	return sum([i*j for (i, j) in zip(u, v)])

def angle_between_vectors(u, v):
	""" Compute the angle between vector u and v """
	dp = dot_product(u, v)
	if dp == 0:
		return 0
	um = magnitude(u)
	vm = magnitude(v)
	return math.acos(dp / (um*vm)) * (180. / math.pi)

def vector_projection(u, v): # projection of u onto v
	scale = dot_product(u, v)/dot_product(v, v)
	proj = vector_mult_s(v, scale)
	return proj

def vector_perpendicular(u, v): # component of u perpendicular to v
	proj = vector_projection(u, v)
	perp = vector_subtract(u, proj)
	return perp

def vector_comp(u, v): # combined projection and perpendicular
	v0 = dot_product(u, v)/magnitude(v)
	scale = dot_product(u, v)/dot_product(v, v)
	proj = vector_mult_s(v, scale)
	perp = vector_subtract(u, proj)
	return {'proj': proj, 'perp': perp, 'v0': v0}

## Define more complex functions that use the vector math functions
# for optimization of the Brachistochrone problem

def velocity_error(vessel, tgt_pos, ref, diag=False):
	now = sc.ut
	vessel_pos = vessel.position(ref)
	vessel_vel = vessel.velocity(ref)
	
	# Get the direction vector pointing towards the target (from vessel to target)
	# dir_vect = tgt_pos - vessel_pos
	# NOTE: Assumes tgt_pos is in reference_frame "ref"
	dir_vect = vector_subtract(tgt_pos, vessel_pos)
	if diag:
		print('Target position: ', tgt_pos)
		print('Vessel position: ', vessel_pos)
		print('Direction Vector:', dir_vect)
		print('Vessel Velocity: ', vessel_vel)
	
	comp = vector_comp(vessel_vel, dir_vect)
	v0 = comp['v0']
	vproj = comp['proj']
	vperp = comp['perp']
	
	if diag:
		print('Initial "helpful" Velcoity:', v0)
		print('"Helpful" velocity vector: ', vproj, magnitude(vproj))
		print('"Harmful" velocity vector: ', vperp, magnitude(vperp))
	
	# Get component of vessel velocity in the direction of the target
	# v0 = dot_product(dir_vect,vessel_vel)/magnitude(dir_vect)
	# A = vessel_vel
	# B = dir_vect
	# scale = dot_product(A, B)/dot_product(B, B)
	# scale = dot_product(A, B)/magnitude(B)**2.0
	# vproj = vector_mult_s(B, scale)
	# vproj = vector_projection(vessel_vel,dir_vect)
	
	# Get components of vessel velocity perpendicular to the direction of the target
	# vperp = vector_subtract(vessel_vel, vproj)
	
	# Get the distance to the target
	dist = magnitude(dir_vect)
	
	info = {"valid_time": now, "dir_vect": dir_vect, "distance": dist, "v0": v0,
	"vperp": vperp, "vproj": vproj}
	
	if diag:
		print(info)
	
	return info
	
def find_optimal_sln(vessel, target, ref, initial_guess=None):
	# Target acceleration in g's - governs thrust calculations
	acc_tgt_g = 10
	acc_tgt = acc_tgt_g*kerbin.surface_gravity
	
	# Compute flip time?
	flip_time = 5.0
	
	print('*** Finding optimal solution for travel to %s ***' % target.name)
	# Get target and vessel current position in sun_rf
	tgt_pos = target.position(ref)
	vessel_pos = vessel.position(ref)
	vessel_vel = vessel.velocity(ref)
	
	# Assume large position error to kick off interration
	position_err = (1e6,1e6,1e6)
	threshold_err = 1e4
	if initial_guess is None:
		tgt_final_pos = tgt_pos
	else:
		tgt_final_pos = initial_guess
	
	passnum = 0
	
	while magnitude(position_err) > threshold_err:
		passnum += 1
		
		# Get current velocity and distance info
		info = velocity_error(vessel, tgt_final_pos, ref)
		now = info['valid_time']
		dir_vect = info['dir_vect']
		dist = info['distance']
		v0 = info['v0']
		vproj = info['vproj']
		vperp = info['vperp']
		
		# print('distance = ', dist)
		
		d3 = (0.5*v0**2.0)/acc_tgt
		dt3 = v0/acc_tgt
		
		if v0 < 0.0:
			d2 = (dist-d3)/2.0
			dt2 = math.sqrt(2.0*d2/acc_tgt)
			d1 = d2 +  d3
			dt1 = dt2 + dt3
		else:
			d1 = (dist-d3)/2.0
			dt1 = math.sqrt(2*d1/acc_tgt)
			d2 = d1 +  d3
			dt2 = dt1 + dt3
		
		# Get the time needed to travel that distance via torch ship burns
		travelTime = dt1 + dt2 + flip_time
		# print('travelTime = ', str(datetime.timedelta(seconds=travelTime)))
		
		# General equation for distance traveled under constant acceleration
		# d = v0*t + 0.5a*t**2
		
		# distance traveled under first part of journey (t0 - t1)
		# t0 = now or start of journey
		# t1 = point of flip
		# t2 = time at end of journey
		# dt1 = t1 - t0 = t1
		# dt2 = t2 - t1
		
		# v0 = initial velocity
		# v1 = velocity at point of flip
		# v2 = final velocity at end of deceleration burn, generally 0?
		# d1 = v0*dt1 + 0.5a*dt1**2
		
		# velocity at point of flip (t1)
		# v1 = v0 + a*t1
		# d2 = v1*dt2 - 0.5*a*dt2**2
		# d1 + d2 = dist -> d2 = dist - d1
		# dist - d1 = v1*dt2 - 0.5*a*dt2**2
		# d1 = dist - v1*dt2 + 0.5*a*dt2**2
		# v2 = v1 -a*dt2 = 0
		# v1 = a*dt2
		
		# travelTime = 2.0*math.sqrt(dist/acc_tgt)
		
		# Get the target position after the burn
		tgt_new_pos = target.orbit.position_at(now + travelTime, ref)
		# print('tgt_new_pos = ', tgt_new_pos)
		
		# get the position error on this iteration
		position_err = vector_subtract(tgt_new_pos, tgt_final_pos)
		# print('%d position_err = %.1f' % (passnum, magnitude(position_err)))
		
		# Update tgt_final_pos
		tgt_final_pos = tgt_new_pos
		# print('tgt_final_pos = ', tgt_final_pos)
				
	sln ={"valid_time": now, "dir_vect": dir_vect, "distance": dist,
	"travel_time": travelTime, "pos_err": position_err, "passnum": passnum,
	"tgt_pos": tgt_final_pos, "dt1": dt1, "dt2": dt2, "d1": d1, "d2": d2,
	"vperp": vperp, "vproj": vproj, "v0": v0, "found":True}
		
	# To really make this optimal we should work out when to depart from the current orbit
	# so that our orbital velocity vector at departure is best aligned with dir_vect, then
	# re-optimize at that time just to be sure we're still within tolerance
	
	# Update the trip info
	# update_rTime_txt(0.0, travelTime, True)
	# update_rDist_txt(dStart, dStart, True)
	
	return sln

def make_node(vessel, t, dV, ref):
	# add_node(ut[, prograde = 0.0][, normal = 0.0][, radial = 0.0])
	# Creates a maneuver node at the given universal time, and returns a Node object that can be used to modify it. Optionally sets the magnitude of the delta-v for the maneuver node in the prograde, normal and radial directions.

	# Parameters:	(in Vessel Orbital Reference Frame)
	# 	ut (double) – Universal time of the maneuver node.
	# 	prograde (float) – Delta-v in the prograde direction.
	# 	normal (float) – Delta-v in the normal direction.
	# 	radial (float) – Delta-v in the radial direction.
	# Return type:	Node
	# Game Scenes:	Flight
	
	# Get the current vessel orbital reference frame
	ves_ref = vessel.orbital_reference_frame
	
	# Get the current vessel position in the same frame as dV
	ves_pos = vessel.position(ref)
	
	# Convert dV from reference_frame ref to the vessel orbital reference frame
	(prograde, normal, radial) = sc.transform_velocity(ves_pos, dV, ref, ves_ref)
	
	# Creaste the node
	node = vessel.control.add_node(t, prograde, normal, radial)
	
	return node

## Connect to tne server and setup the UI
#This script assumes the vessel is in orbit and the target is set.
conn = krpc.connect(name='TorchShipControl')

# Get Space Center and MechJeb objects
sc = conn.space_center
mj = conn.mech_jeb

# Get the canvas we'll place the panel on
canvas = conn.ui.stock_canvas

# Get the size of the game window in pixels
screen_size = canvas.rect_transform.size

# Add a panel to contain the UI elements
panel = canvas.add_panel()

# Position the panel on the left of the screen
rect = panel.rect_transform
num_rows = 5
font_size = 14
margin = 15
x_size = 300
y_size = num_rows*20 + 20
rect.size = (x_size, y_size)
# Position is measured from the center of the screen_size
# X range is screen_size[0]*(-1/2, 1/2)
# Y range is screen_size[1]*(-1/2, 1/2)
# Centered horizontally at 200 pixels from left edge, and centered vertically on screen
x_pos = (x_size - 100) - screen_size[0]/2
y_pos = 0
rect.position = (x_pos, y_pos)

row = 1
y_pos = (num_rows/2.0 -row)*20

w1 = 70
h = 30 # 25?
x_pos = margin + w1/2.0 - x_size/2.0

# Add a button to control the Torch Pilot action (run/stop)
actButton = panel.add_button("Run")
actButton.rect_transform.position = (x_pos - 5, y_pos + 10)
actButton.rect_transform.size = (w1, h - 5)

w2 = 50
x_pos += w1/2.0 + w2/2.0 + 10

# Add an input field for the target
target_label = panel.add_text('Target:')
target_label.rect_transform.position = (x_pos, y_pos + 5) # was y_pos + 5
target_label.rect_transform.size = (w2, h)
target_label.color = (1, 1, 1)
target_label.size = font_size+1

w3 = 130
x_pos += w2/2.0 + w3/2.0 + 10

target_txt = panel.add_text('Unavailable')
target_txt.rect_transform.position = (x_pos, y_pos + 5) # was y_pos + 5
target_txt.rect_transform.size = (w3, h)
target_txt.color = (1, 1, 1)
target_txt.size = font_size+1

# Unable to get target_body in 'space_center', 'tracking_station', 'editor_sph', 'editor_vab', 
if conn.krpc.current_game_scene.name == "flight":
	if sc.target_body is None:
		target_txt.content = 'Unset'
	else:
		target_txt.content = sc.target_body.name

#new_target = panel.add_input_field()
#new_target.rect_transform.position = (70, 30)
#new_target.rect_transform.size = (140, 25)
# Unable to get target_body in 'space_center', 'tracking_station', 'editor_sph', 'editor_vab', 
#if conn.krpc.current_game_scene.name = "flight":
#	new_target.value = sc.target_body.name
#else:
#	new_target.value = 'unavailable'

# Set up a stream to monitor the throttle actButton
actButton_clicked = conn.add_stream(getattr, actButton, 'clicked')

row = 2
y_pos = (num_rows/2.0 -row)*20

w = 270
# h = 30
x_pos = margin + w/2.0 - x_size/2.0

# Add some text displaying the total engine thrust
rTime_txt = panel.add_text("Remaining Time: 00:00:00 (100.0%)")
rTime_txt.rect_transform.size = (w, h)
rTime_txt.rect_transform.position = (x_pos, y_pos)
rTime_txt.color = (1, 1, 1)
rTime_txt.size = font_size
# rTime_txt.alignment

row = 3
y_pos = (num_rows/2.0 -row)*20

w = 270
# h = 30
x_pos = margin + w/2.0 - x_size/2.0

rDist_txt = panel.add_text("Remaining Distance: 0 km (100.0%)")
rDist_txt.rect_transform.size = (w, h)
rDist_txt.rect_transform.position = (x_pos, y_pos)
rDist_txt.color = (1, 1, 1)
rDist_txt.size = font_size

row = 4
y_pos = (num_rows/2.0 -row)*20

w1 = 35
# h = 30
x_pos = margin + w1/2.0 - x_size/2.0

hErr_lable = panel.add_text("H Err:")
hErr_lable.rect_transform.size = (w1, h)
hErr_lable.rect_transform.position = (x_pos, y_pos)
hErr_lable.color = (1, 1, 1)
hErr_lable.size = font_size

w2 = 105
x_pos += w1/2.0 + w2/2.0 + 10

hErr_txt = panel.add_text("-000.000" u"\N{DEGREE SIGN}")
hErr_txt.rect_transform.size = (w2, h)
hErr_txt.rect_transform.position = (x_pos, y_pos)
hErr_txt.color = (1, 1, 1)
hErr_txt.size = font_size

w3 = 35
x_pos += w2/2.0 + w3/2.0 + 10

vErr_lable = panel.add_text("V Err:")
vErr_lable.rect_transform.size = (w3, h)
vErr_lable.rect_transform.position = (x_pos, y_pos)
vErr_lable.color = (1, 1, 1)
vErr_lable.size = font_size

w4 = 75
x_pos += w3/2.0 + w4/2.0 + 10

vErr_txt = panel.add_text("-000.0 m/s")
vErr_txt.rect_transform.size = (w4, h)
vErr_txt.rect_transform.position = (x_pos, y_pos)
vErr_txt.color = (1, 1, 1)
vErr_txt.size = font_size

row = 5
y_pos = (num_rows/2.0 -row)*20

w1 = 45
x_pos = margin + w1/2.0 - x_size/2.0

status_lable = panel.add_text("Status:")
status_lable.rect_transform.size = (w1, h)
status_lable.rect_transform.position = (x_pos, y_pos)
status_lable.color = (1, 1, 1)
status_lable.size = font_size

w2 = 95
x_pos += w1/2.0 + w2/2.0 + 10

status_txt = panel.add_text("Idle")
status_txt.rect_transform.size = (w2, h)
status_txt.rect_transform.position = (x_pos, y_pos)
status_txt.color = (1, 1, 1)
status_txt.size = font_size
status_txt.content = 'Nulling Velocity Error'

w3 = 35
x_pos += w2/2.0 + w3/2.0 + 10

pErr_lable = panel.add_text("P Err:")
pErr_lable.rect_transform.size = (w3, h)
pErr_lable.rect_transform.position = (x_pos, y_pos)
pErr_lable.color = (1, 1, 1)
pErr_lable.size = font_size

w4 = 65
x_pos += w3/2.0 + w4/2.0 + 10

pErr_txt = panel.add_text("-000.000" u"\N{DEGREE SIGN}")
pErr_txt.rect_transform.size = (w4, h)
pErr_txt.rect_transform.position = (x_pos, y_pos)
pErr_txt.color = (1, 1, 1)
pErr_txt.size = font_size

if conn.krpc.current_game_scene.name == "space_center":
	status_txt.content = "unavailable"

# rTime_txt.alignement
# rTime_txt.available_fonts
# rTime_txt.color
# rTime_txt.content
# rTime_txt.font = 'Ariel'
# rTime_txt.line_spacing = 1.0
# rTime_txt.rect_transform.anchor
# rTime_txt.rect_transform.anchor_max = (0.5, 0.5)
# rTime_txt.rect_transform.anchor_min = (05., 0.5)
# rTime_txt.rect_transform.local_position = (-60.0, 5.0, 0.0)
# rTime_txt.rect_transform.lower_left = (-140.0, -10.0)
# rTime_txt.rect_transform.pivot = (0.5, 0.5)
# rTime_txt.rect_transform.position = (-60.0, 5.0)
# rTime_txt.rect_transform.rotation = (0.0, 0.0, 0.0, 1.0)
# rTime_txt.rect_transform.scale = (1.0, 1.0, 1.0)
# rTime_txt.rect_transform.size = (160.0, 30.0)


# Function to update the elapsed time in the journey
def update_rTime_txt(dt, tt, diag=False):
	frac = (tt-dt)/tt
	et = str(datetime.timedelta(seconds=round(tt-dt)))
	s = 'Remaining Time: %s (%.3f%%)' % (et, 100.0*frac)
	# rTime_txt.content = 'Elapsed Time: %s (%.3f%%)' % (et, 100.0*frac)
	rTime_txt.content = s
	if diag:
		print(s)

# Function to update the remaining distance to the target
def update_rDist_txt(dist, td, diag=False):
	frac = dist/td
	if dist > 1e9:
		dist /= 1e9
		units = 'Gm'
	elif dist > 1e6:
		dist /= 1e6
		units = 'Mm'
	elif dist > 1e3:
		dist /= 1e3
		units = 'km'
	else:
		units = 'm'
	s = 'Remaining Distance: %.1f%s (%.3f%%)' % (dist, units, 100.0*frac)
	# rDist_txt.content = 'Distance: %.1f%s (%.3f%%)' % (dist, units, 100.0*frac)
	rDist_txt.content = s
	if diag:
		print(s)

# Function to update the remaining distance to the target
def update_hErr_txt(h_err, diag=False):
	s = '%.3f' % h_err + u'\N{DEGREE SIGN}'
	# rDist_txt.content = 'Distance: %.1f%s (%.3f%%)' % (dist, units, 100.0*frac)
	hErr_txt.content = s
	if diag:
		print("Heading Err:", s)

# Function to update the remaining distance to the target
def update_vErr_txt(v_err, diag=False):
	if v_err > 1e9:
		v_err /= 1e9
		units = 'Gm/s'
	elif v_err > 1e6:
		v_err /= 1e6
		units = 'Mm/s'
	elif v_err > 1e3:
		v_err /= 1e3
		units = 'km/s'
	else:
		units = 'm/s'
	s = '%.1f %s' % (v_err, units)
	vErr_txt.content = s
	if diag:
		print("Magnitude of Velocty Err:", s)

# Function to update the pointing error to the target
def update_pErr_txt(p_err, diag=False):
	s = '%.3f' % p_err + u'\N{DEGREE SIGN}'
	pErr_txt.content = s
	if diag:
		print("Pointing Err:", s)

# Function to update the current status
def update_status_txt(state, diag=False):
	if state == 1:
		status = 'Idle'
		# status_txt.color = [0, 0, 1]
	elif state == 2:
		status = 'Nulling Velocity Error'
		# status_txt.color = [0, 1, 0]
	elif state == 3:
		status = 'Pointing to %s' % sc.target_body.name
		# status_txt.color = [0, 1, 0]
	elif state == 4:
		status = 'Burn, Baby! Burn!'
		# status_txt.color = [0, 0, 1]
	elif state == 5:
		status = 'Flippin!'
		# status_txt.color = [0, 1, 0]
	elif state == 6:
		status = 'I Said Stop!'
		# status_txt.color = [1, 0, 0]
	elif state == 7:
		status = 'Set Inclination'
		# status_txt.color = [0, 1, 0]
	elif state == 8:
		status = 'Set Periapsis'
		# status_txt.color = [0, 1, 0]
	elif state == 9:
		status = 'Circularize!'
		# status_txt.color = [0, 1, 0]
	elif state == 10:
		status = 'Done!'
		# status_txt.color = [0, 1, 0]
	else:
		status = 'Unset'
	
	status_txt.content = status
	if diag:
		print('Status: %s' % status)

executor = mj.node_executor
def execute_nodes():
	print("Executing maneuver nodes")
	executor.execute_all_nodes()
	
	with conn.stream(getattr, executor, "enabled") as enabled:
		# We don't need a high throughput rate, 1 second is more than enough
		enabled.rate = 1
		with enabled.condition:
			while enabled():
				enabled.wait()

# Get the Sun object from the sc.bodies dictionary
sun = sc.bodies['Sun']
sun_rf = sun.non_rotating_reference_frame

# Get other celestial bodies we may need
moho = sc.bodies['Moho']
eve = sc.bodies['Eve']
gilly = sc.bodies['Gilly']
kerbin = sc.bodies['Kerbin']
mun = sc.bodies['Mun']
minmus = sc.bodies['Minmus']
duna = sc.bodies['Duna']
ike = sc.bodies['Ike']
dres = sc.bodies['Dres']
jool = sc.bodies['Jool']
laythe = sc.bodies['Laythe']
vall = sc.bodies['Vall']
tylo = sc.bodies['Tylo']
bop = sc.bodies['Bop']
pol = sc.bodies['Pol']
eeloo = sc.bodies['Eeloo']

node = None
t_line = None
initial_path_line_r1 = None
initial_path_line_r2 = None
final_path_line_r1 = None
final_path_line_r2 = None
vproj_line_r1 = None
vperp_line_r1 = None
dvec_tt = None
dvec_vt = None
dvec_ve = None
dvec_vv = None

# Initialize the Torch Pilot state
tp_state = 1
update_status_txt(tp_state)

ref = celestial_body_nr_rf
# ref = sun_rf
ref2 = orbit_rf
# ref2 = vessel_rf

red = (1,0,0)
green = (0,1,0)
blue = (0,0,1)
yellow = (1,1,0)
magenta = (1,0,1)
cyan = (0,1,1)
white = (1,1,1)

# Endless outer loop
while True:
	# Process actButton press
	if actButton_clicked():
		conn.ui.message('TorchPilot only runs in Flight Mode!', color=(1,0,0))
	
	sass = None
	vessel = None
	target = None
	# If we're in Flight mode, then enter endless loop for that mode
	while conn.krpc.current_game_scene.name == "flight":
		# Do first time stuff - doesn't change with vessel changes
		if sass is None:
			# Get the SmartASS module object for MechJeb (only works in Flight mode)
			sass = mj.smart_ass
			
		# If the target has changed or is not yet defined
		if target != sc.target_body:
			# Record the "new" target object
			target = sc.target_body
			
			# Report it on the GUI
			if target is None:
				target_txt.content = 'Unset'
			else:
				target_txt.content = sc.target_body.name
				
				# Require a new optimal solution
				need_optimal_sln = True
		
		# If the active vessel has changed or is not yet defined 
		if vessel != sc.active_vessel:
			# Record the active vessel object
			vessel = sc.active_vessel
			surf_vel_rf = vessel.surface_velocity_reference_frame
			orbit_rf = vessel.orbital_reference_frame
			celestial_body_rf = vessel.orbit.body.reference_frame
			celestial_body_nr_rf = vessel.orbit.body.non_rotating_reference_frame
			vessel_rf = vessel.reference_frame


			# Get the vessel auto_pilot for later use
			ap = vessel.auto_pilot
			ap.reference_frame = sun_rf

			# Require a new optimal solution
			if target is None:
				need_optimal_sln = False
			else:
				need_optimal_sln = True

		vessel_pos = vessel.position(sun_rf)

		# If we need a new optimal solution...
		if need_optimal_sln:
			sln = find_optimal_sln(vessel, target, sun_rf)
			if sln["found"]:
				print('Solution found after %d itterations! (pos_err = %.1f)' % (sln['passnum'], magnitude(sln['pos_err'])))
				print(sln)
				need_optimal_sln = False
				print('Travel Time = ', str(datetime.timedelta(seconds=sln['travel_time'])))
				print('Travel Distance = %.1f' % sln['distance'])
				h_err = angle_between_vectors(sln['dir_vect'], vessel.velocity(sun_rf))
				update_hErr_txt(h_err, True)
				update_vErr_txt(magnitude(sln['vperp']), True)
				p_err = angle_between_vectors(sln['dir_vect'], vessel.direction(sun_rf))
				update_pErr_txt(p_err, True)
				
				if not(initial_path_line is None):
					initial_path_line.remove()
				tgt_pos = target.position(sun_rf)
				print('Drawing line from vessel to target')
				initial_path_line = conn.drawing.add_line(vessel_pos, tgt_pos, sun_rf)
				initial_path_line.color = (0,0,1)
				print('\tthickness = ', initial_path_line.thickness)
				print('\tcolor = ', initial_path_line.color)
				print('\tvisible = ', initial_path_line.visible)
				
				if not(final_path_line is None):
					final_path_line.remove()
				tgt_final_pos = sln['tgt_pos']
				print('Drawing line from vessel to virtual target')
				final_path_line = conn.drawing.add_line(vessel_pos, tgt_final_pos, sun_rf)
				final_path_line.color = (0,0,1)
				print('\tthickness = ', final_path_line.thickness)
				
				vproj = sln['vproj']
				vperp = sln['vperp']
				vproj_end = vector_add(vessel_pos, vproj)
				vperp_end = vector_add(vessel_pos, vperp)
				
				if not(vproj_line is None):
					vproj_line.remove()
				print('Drawing vproj line from vessel')
				vproj_line = conn.drawing.add_line(vessel_pos, vproj_end, sun_rf)
				vproj_line.color = (0,1,0)
				vproj_line.thickness *= 2.0
				
				if not(vperp_line is None):
					vperp_line.remove()
				print('Drawing vproj line from vessel')
				vperp_line = conn.drawing.add_line(vessel_pos, vperp_end, sun_rf)
				vperp_line.color = (1,0,0)
				vperp_line.thickness *= 2.0
				
		if not(need_optimal_sln):
			# Get current velocity and distance info
			info = velocity_error(vessel, sln['tgt_pos'], sun_rf)
			now = info['valid_time']
			dir_vect = info['dir_vect']
			dist = info['distance']
			v0 = info['v0']
			vproj = info['vproj']
			vperp = info['vperp']
			
			td = sln['distance']
			update_rDist_txt(dist, td)
			
			tt = sln['travel_time']
			if actButton.text.content == 'Run':
				dt = 0.0
			else:
				dt = sc.ut - sln["valid_time"]
			update_rTime_txt(dt, tt)

			h_err = angle_between_vectors(dir_vect, vessel.velocity(sun_rf))
			update_hErr_txt(h_err)
			
			update_vErr_txt(magnitude(vperp))

			if node is None:
				p_err = angle_between_vectors(sln['dir_vect'], vessel.direction(sun_rf))
				# print("pointing Error (No Node): ", p_err)
			else:
				p_err = angle_between_vectors(sln['dir_vect'], node.burn_vector(sun_rf))
				# print("pointing Error (With Node): ", p_err)
			update_pErr_txt(p_err)
				
		# Process actButton press
		if actButton_clicked():
			actButton.clicked = False
			if actButton.text.content == 'Run':
				if need_optimal_sln:
					sln = find_optimal_sln(vessel, target, sun_rf)
					if sln["found"]:
						print('Solution found after %d itterations! (pos_err = %.1f)' % (sln['passnum'], magnitude(sln['pos_err'])))
						# print('Solution found!')
						print(sln)
						need_optimal_sln = False
				elif sc.ut - sln["valid_time"] > 100:
					need_optimal_sln = True
					sln = find_optimal_sln(vessel, target, sun_rf, sln['tgt_pos'])
					if sln["found"]:
						print('New solution found after %d itterations! (pos_err = %.1f)' % (sln['passnum'], magnitude(sln['pos_err'])))
						print(sln)
						need_optimal_sln = False
				print("Commencing Run...")
				tp_state = 2
				node = None
				actButton.text.content = 'Stop'
			else:
				print("Aborting Run...")
				# Kill Thrust!
				tp_state = 1
				node = None
				actButton.text.content = 'Run'
				
		if (tp_state > 1) and not need_optimal_sln:
			# If the tp_state is > 1, then an optimal solution must exist and we can get
			# the current velocity error, etc.
			info = velocity_error(vessel, sln['tgt_pos'], sun_rf)
			vperp = info['vperp']
			dist = info['distance']
			
		if tp_state == 2:
			# State: Nulling Velocity Error
			if (magnitude(vperp) > 1) and node is None:
				dV = vector_mult_s(vperp, -1.0)
				node = make_node(vessel, sc.ut, dV, sun_rf)
				sass.autopilot_mode = sass.autopilot_mode.node
				sass.update(False)
				execute_nodes()
				# sass.autopilot_mode = sass.autopilot_mode.off
				# tp_state = 3
			
			# Command vessel to point away from vperp
			# if heading error is small, then apply thrust
			# if vperp is small cut thrust and advance state
		elif tp_state == 3:
			print('Entering State 3')
			# status = 'Pointing to %s' % sc.target_body.name
			# Command vessel to point at dir_vect
			# if heading error is small, then advance state
		elif tp_state == 4:
			print('Entering State 4')
			# status = 'Burn, Baby! Burn!'
			# if heading error is small and dt < t1, then set thrust on
		elif tp_state == 5:
			print('Entering State 5')
			# status = 'Flippin!'
			# Command vessel to point at -dir_vect
			# if heading error is small, then advance state
		elif tp_state == 6:
			print('Entering State 6')
			# status = 'I Said Stop!'
			# if heading error is small and dt < t2, then set thrust on
			# if dt >= t2, then set thrust off and advance state
		elif tp_state == 7:
			print('Entering State 7')
			# status = 'Set Inclination'
			# use MechJeb to set inclination to zero-ish
			# if done, then advance state
		elif tp_state == 8:
			print('Entering State 8')
			# status = 'Set Periapsis'
			# use MechJeb st set the periapsis to something reasonable
			# if done, then advance state
		elif tp_state == 9:
			print('Entering State 9')
			# status = 'Circularize!'
			# use MechJeb to circularize at periapsis
			# If done, then advance state
		elif tp_state > 9:
			print("Done!")
			tp_state = 1
			
		update_status_txt(tp_state)
		
		# If we're not already orbiting the target_body, find the optimal solution
		# Find optimal solution
		
		# Wait for optimal departure time
		# Burn to remove launch errors
		# Accelerate along optimum trajectory toward target
		# Decelerate to arrive at futur target location
		# Match velocity with target
		# Use MJ to get into the fina orbit
		
		# wait 0.1 seconds before next pass through the inner loop
		time.sleep(0.1)
		
	# Wait 0.1 seconds before next pass through the outer loop	
	time.sleep(0.1)

print("Done!")
conn.close()
