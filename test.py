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
		info = velocity_error(vessel, tgt_final_pos, ref, True)
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
# mj = conn.mech_jeb

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

# sass = None
vessel = None
target = None

vessel = conn.space_center.active_vessel

surf_vel_rf = vessel.surface_velocity_reference_frame
orbit_rf = vessel.orbital_reference_frame
celestial_body_rf = vessel.orbit.body.reference_frame
celestial_body_nr_rf = vessel.orbit.body.non_rotating_reference_frame
vessel_rf = vessel.reference_frame

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

dvec_x = conn.drawing.add_direction((1, 0, 0), orbit_rf)
dvec_y = conn.drawing.add_direction((0, 1, 0), orbit_rf)
# dvec_z = conn.drawing.add_direction((0, 0, 1), orbit_rf)
dvec_x.color = red
dvec_y.color = green
# dvec_z.color = blue

sln_time_tol = 60.0
first_pass = True

# Endless outer loop
while True:
	if target != sc.target_body:
		# Record the "new" target object
		print('Updating Target')
		target = sc.target_body
		need_optimal_sln = True

	while target == sc.target_body:
		tgt_pos = target.position(ref)
		vessel_pos = vessel.position(ref)
		# if t_line is None:
		# 	t_line = conn.drawing.add_line(vessel_pos, tgt_pos, ref)
		# 	t_line.color = blue
		# 	# t_line.thickness = 0.25
		# else:
		# 	t_line.start = vessel_pos
		# 	t_line.end = tgt_pos
		
		if need_optimal_sln:
			sln = find_optimal_sln(vessel, target, ref)
			if sln["found"]:
				print('')
				print('Solution found after %d itterations! (pos_err = %.1f)' % (sln['passnum'], magnitude(sln['pos_err'])))
				print(sln)
				print('Travel Time = ', str(datetime.timedelta(seconds=sln['travel_time'])))
				print('Travel Distance = %.1f' % sln['distance'])
				h_err = angle_between_vectors(sln['dir_vect'], vessel.velocity(ref))
				print('h_err = ', h_err)
				sgn = 1
				if h_err > 90:
					sgn = -1
				print('vproj = ', sln['vproj'], sgn*magnitude(sln['vproj']))
				print('vperp = ', sln['vperp'], magnitude(sln['vperp']))
				need_optimal_sln = False
				first_pass = True
			
		## Record positions at time now
		# Current target position
		tgt_pos_r1 = target.position(ref)
		tgt_pos_r2 = target.position(ref2)
		# Current vessel position
		vessel_pos_r1 = vessel.position(ref)
		vessel_pos_r2 = vessel.position(ref2)
		# Virtual target position (target will be here when we arrive)
		vtgt_pos_r1 = sln['tgt_pos'] # propagate this?
		vtgt_pos_r2 = sc.transform_position(vtgt_pos_r1, ref, ref2)
		
		## Record vectors at time now
		# Current target direction
		tgt_dir_r1 = vector_subtract(tgt_pos_r1, vessel_pos_r1)
		tgt_dir_r2 = vector_subtract(tgt_pos_r2, vessel_pos_r2)
		# Direction from vessel to virtual target
		vtgt_dir_r1 = sln['dir_vect'] # propagate this?
		vtgt_dir_r2 = sc.transform_direction(vtgt_dir_r1, ref, ref2)
		# Direction of vessel current velocity
		vessel_vel_r1 = vessel.velocity(ref)
		# vessel_vel_r2 = vessel.velocity(ref2)
		vessel_vel_r2 = sc.transform_direction(vessel_vel_r1, ref, ref2)
		
		# Direction of vessel body orientation
		vessel_body_r1 = vessel.direction(ref)
		vessel_body_r2 = vessel.direction(ref2)
		# Direction of vessel velocity error
		vperp_r1 = sln['vperp']
		# vperp_r2 = sc.transform_velocity(vessel_pos_r1, vperp_r1, ref, ref2)
		vperp_r2 = sc.transform_direction(vperp_r1, ref, ref2)
		# Direction of vessel velocity error
		vproj_r1 = sln['vproj']
		# vperp_r2 = sc.transform_velocity(vessel_pos_r1, vproj_r1, ref, ref2)
		vproj_r2 = sc.transform_direction(vproj_r1, ref, ref2)
				
		if first_pass:
			# Display heading error (angular diff between vtgt_dir and vessel velocity vector)
			h_err_r1 = angle_between_vectors(vtgt_dir_r1, vessel_vel_r1)
			h_err_r2 = angle_between_vectors(vtgt_dir_r2, vessel_vel_r2)
			# update_hErr_txt(h_err, True)
			s1 = '%.3f' % h_err_r1 + u'\N{DEGREE SIGN}'
			# s2 = '%.3f' % h_err_r2 + u'\N{DEGREE SIGN}'
			print("Heading Error:", s1)
			sgn = 1
			if h_err_r1 > 90:
				sgn = -1

			# Display pointing error (angular diff between vtgt_dir and vessel body orientation)
			p_err_r1 = angle_between_vectors(vtgt_dir_r1, vessel_body_r1)
			p_err_r2 = angle_between_vectors(vtgt_dir_r2, vessel_body_r2)
			# update_pErr_txt(p_err, True)
			s1 = '%.3f' % p_err_r1 + u'\N{DEGREE SIGN}'
			s2 = '%.3f' % p_err_r2 + u'\N{DEGREE SIGN}'
			print("Pointing Error:", s1, s2)
			
			print("ref:")
			print("tgt_dir_r1:   ", tgt_dir_r1, magnitude(tgt_dir_r1))
			print("vessel_vel_r1:", vessel_vel_r1, magnitude(vessel_vel_r1))
			print("vtgt_dir_r1:  ", vtgt_dir_r1, magnitude(vtgt_dir_r1))
			print("vproj_r1:     ", vproj_r1, sgn*magnitude(vproj_r1))
			print("vperp_r1:     ", vperp_r1, magnitude(vperp_r1))
			print("ref2:")
			print("tgt_dir_r2:   ", tgt_dir_r2, magnitude(tgt_dir_r2))
			print("vessel_vel_r2:", vessel_vel_r2, magnitude(vessel_vel_r2))
			print("vtgt_dir_r2:  ", vtgt_dir_r2, magnitude(vtgt_dir_r2))
			print("vproj_r2:     ", vproj_r2, sgn*magnitude(vproj_r2))
			print("vperp_r2:     ", vperp_r2, magnitude(vperp_r2))
			
		# Draw direction line for Vessel -> Target (white)
		if not(dvec_tt is None):
			dvec_tt.remove() 
		dvec_tt = conn.drawing.add_direction(tgt_dir_r2, ref2)
		dvec_tt.color = white
	
		# Draw direction line for Vessel Velocity (red)
		if not(dvec_vv is None):
			dvec_vv.remove() 
		# dvec_vv = conn.drawing.add_direction(vessel_vel_r1, ref)
		# dvec_vv.color = magenta
		dvec_vv = conn.drawing.add_direction(vessel_vel_r2, ref2)
		dvec_vv.color = magenta
	
		# Draw direction line for Vessel -> Virtual Target (cyan)
		if not(dvec_vt is None):
			dvec_vt.remove() 
		# dvec_vt = conn.drawing.add_direction(vtgt_dir_r1, ref)
		# dvec_vt.color = cyan
		dvec_vt = conn.drawing.add_direction(vtgt_dir_r2, ref2)
		dvec_vt.color = cyan
	
		# update_vErr_txt(magnitude(sln['vperp']), True)
		# Draw direction line for Vessel Velocity Error (yellow)
		if not(dvec_ve is None):
			dvec_ve.remove() 
		# dvec_ve = conn.drawing.add_direction(vperp_r1, ref)
		# dvec_ve.color = yellow
		dvec_ve = conn.drawing.add_direction(vector_mult_s(vperp_r2, -1.0), ref2)
		dvec_ve.color = yellow
			
		# else:
		# 	dvec_tt.start = vessel_pos_r1
		# 	dvec_vv.start = vessel_pos_r1
		# 	dvec_vt.start = vessel_pos_r1
		# 	dvec_ve.start = vessel_pos_r1
			
		# if initial_path_line_r1 is None:
		# 	initial_path_line_r1 = conn.drawing.add_line(vessel_pos_r1, tgt_pos_r1, ref)
		# 	initial_path_line_r1.color = (1,0,0)
		# 	initial_path_line_r1.thickness = 0.5
		# else:
		# 	initial_path_line_r1.start = vessel_pos_r1
		# 	initial_path_line_r1.end = tgt_pos_r1
			
		# if initial_path_line_r2 is None:
		# 	initial_path_line_r2 = conn.drawing.add_line(vessel_pos_r2, tgt_pos_r2, ref2)
		# 	initial_path_line_r2.color = (0,1,1)
		# 	initial_path_line_r2.thickness = 0.5
		# else:
		# 	initial_path_line_r2.start = vessel_pos_r2
		# 	initial_path_line_r2.end = tgt_pos_r2

		
		# if final_path_line_r1 is None:
		# 	final_path_line_r1 = conn.drawing.add_line(vessel_pos_r1, vtgt_pos_r1, ref)
		# 	final_path_line_r1.color = (0,0,1)
		# 	final_path_line_r1.thickness = 0.5
		# else:
		# 	final_path_line_r1.start = vessel_pos_r1
		# 	final_path_line_r1.end = vtgt_pos_r1
		
		# if final_path_line_r2 is None:
		# 	final_path_line_r2 = conn.drawing.add_line(vessel_pos_r2, vtgt_pos_r2, ref2)
		# 	final_path_line_r2.color = (0,1,0)
		# 	final_path_line_r2.thickness = 0.5
		# else:
		# 	final_path_line_r1.start = vessel_pos_r2
		# 	final_path_line_r1.end = vtgt_pos_r2
		
		# print('\tthickness = ', final_path_line.thickness)
		
		# dvt = vector_subtract(tgt_final_pos, vp2)
		# dvec_vt = conn.drawing.add_direction(dvt, ref)
		# dvec_vt.color = (1,1,1)
		# 
		# vproj = sc.transform_position(sln['vproj'], ref, ref2)
		# vperp = sc.transform_position(sln['vperp'], ref, ref2)
		# vproj_end = vector_add(vp2, vproj)
		# vperp_end = vector_add(vp2, vperp)
		# 
		# if not(vproj_line is None):
		# 	vproj_line.remove()
		# print('Drawing vproj line from vessel (%.1f long)' % magnitude(vproj))
		# vproj_line = conn.drawing.add_line(vp2, vproj_end, ref2)
		# vproj_line.color = (0,1,1)
		# vproj_line.thickness *= 2.0
		# 
		# if not(vperp_line is None):
		# 	vperp_line.remove()
		# print('Drawing vproj line from vessel (%.1f long)' % magnitude(vperp))
		# vperp_line = conn.drawing.add_line(vp2, vperp_end, ref2)
		# vperp_line.color = (1,0,1)
		# vperp_line.thickness *= 2.0
     	
     	# Check to see if we should trigger finding a new optimal solution
		if sc.ut > (sln['valid_time'] + sln_time_tol):
			need_optimal_sln = True
		first_pass = False

		# wait 0.1 seconds before next pass through the inner loop
		time.sleep(0.1)
			
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
