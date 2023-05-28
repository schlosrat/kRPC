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
	
def cross_product(u, v): # compute vector cross product of u x v
	return (u[1]*v[2] - u[2]*v[1],
	u[2]*v[0] - u[0]*v[2],
	u[0]*v[1] - u[1]*v[0])

def dot_product(u, v): # Compute scalar dot product of u and v
	return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

# def magnitude(v):
# 	return math.sqrt(dot_product(v, v))

def angle_between_vectors(u, v):
	""" Compute the angle between vector u and v """
	dp = dot_product(u, v)
	if dp == 0:
		return 0
	um = magnitude(u)
	vm = magnitude(v)
	return math.acos(dp / (um*vm)) * (180. / math.pi)

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
x_size = 300
y_size = 100
rect.size = (x_size, y_size)
# Position is measured from the center of the screen_size
# X range is screen_size[0]*(-1/2, 1/2)
# Y range is screen_size[1]*(-1/2, 1/2)
# Centered horizontally at 200 pixels from left edge, and centered vertically on screen
x_pos = 200 - screen_size[0]/2
y_pos = 0
rect.position = (x_pos, y_pos)

# Add a button to control the Torch Pilot
button = panel.add_button("Run")
button.rect_transform.position = (-105, 30)
button.rect_transform.size = (70, 25)

# Add an input field for the target
text0 = panel.add_text('Target:')
text0.rect_transform.position = (-20, 25)
text0.rect_transform.size = (70, 25)
text0.color = (1, 1, 1)
text0.size = 14
text_tgt = panel.add_text('Unavailable')
text_tgt.rect_transform.position = (70, 25)
text_tgt.rect_transform.size = (140, 25)
text_tgt.color = (1, 1, 1)
text_tgt.size = 14
# Unable to get target_body in 'space_center', 'tracking_station', 'editor_sph', 'editor_vab', 
if conn.krpc.current_game_scene.name = "flight":
	text_tgt.value = sc.target_body.name

#new_target = panel.add_input_field()
#new_target.rect_transform.position = (70, 30)
#new_target.rect_transform.size = (140, 25)
# Unable to get target_body in 'space_center', 'tracking_station', 'editor_sph', 'editor_vab', 
#if conn.krpc.current_game_scene.name = "flight":
#	new_target.value = sc.target_body.name
#else:
#	new_target.value = 'unavailable'

# Set up a stream to monitor the throttle button
button_clicked = conn.add_stream(getattr, button, 'clicked')

# Add some text displaying the total engine thrust
text1 = panel.add_text("Elapsed Time: 00:00:00 (100.0%)")
text1.rect_transform.size = (270, 30)
text1.rect_transform.position = (0, 0)
text1.color = (1, 1, 1)
text1.size = 14
text1.alignment
text2 = panel.add_text("Remaining Distance: 0 km (100.0%)")
text2.rect_transform.size = (270, 30)
text2.rect_transform.position = (0, -20)
text2.color = (1, 1, 1)
text2.size = 14
text3 = panel.add_text("Status:")
text3.rect_transform.size = (70, 30)
text3.rect_transform.position = (-100, -40)
text3.color = (1, 1, 1)
text4 = panel.add_text("Idle")
text4.rect_transform.size = (220, 30)
text4.rect_transform.position = (25, -40)
text4.color = (1, 1, 1)

if conn.krpc.current_game_scene.name == "space_center":
	text4.content = "unavailable"

# text1.alignement
# text1.available_fonts
# text1.color
# text1.content
# text1.font = 'Ariel'
# text1.line_spacing = 1.0
# text1.rect_transform.anchor
# text1.rect_transform.anchor_max = (0.5, 0.5)
# text1.rect_transform.anchor_min = (05., 0.5)
# text1.rect_transform.local_position = (-60.0, 5.0, 0.0)
# text1.rect_transform.lower_left = (-140.0, -10.0)
# text1.rect_transform.pivot = (0.5, 0.5)
# text1.rect_transform.position = (-60.0, 5.0)
# text1.rect_transform.rotation = (0.0, 0.0, 0.0, 1.0)
# text1.rect_transform.scale = (1.0, 1.0, 1.0)
# text1.rect_transform.size = (160.0, 30.0)


# Function to update the elapsed time in the journey
def update_text1(dt, tt, diag=False):
	frac = dt/tt
	et = str(datetime.timedelta(seconds=dt))
	s = 'Elapsed Time: %s (%.3f%%)' % (et, 100.0*frac)
	# text1.content = 'Elapsed Time: %s (%.3f%%)' % (et, 100.0*frac)
	text1.content = s
	if diag:
		print(s)

# Function to update the remaining distance to the target
def update_text2(dist, td, diag=False):
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
	# text2.content = 'Distance: %.1f%s (%.3f%%)' % (dist, units, 100.0*frac)
	text2.content = s
	if diag:
		print(s)

# Function to update the current status
def update_text4(state, diag=False):
	if state == 1:
		sataus = 'Idle'
		text4.color = [0, 0, 1]
	elif state == 2:
		status = 'Pointing to %s' % sc.target_body.name
		text4.color = [0, 1, 0]
	elif state == 3:
		status = 'Burn, Baby! Burn!'
		text4.color = [0, 0, 1]
	elif state == 4:
		status = 'Flippin!'
		text4.color = [0, 1, 0]
	elif state == 5:
		status = 'I Said Stop!'
		text4.color = [1, 0, 0]
	elif state == 6:
		status = 'Set Inclination'
		text4.color = [0, 1, 0]
	elif state == 7:
		status = 'Set Periapsis'
		text4.color = [0, 1, 0]
	elif state == 8:
		status = 'Circularize!'
		text4.color = [0, 1, 0]

	text4.content = status
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

# Target acceleration in g's - governs thrust calculations
acc_tgt_g = 10
acc_tgt = acc_tgt_g*kerbin.surface_gravity

# Initialize the Torch Pilot state
tp_state = 1
update_text4(tp_state)

# Compute flip time?
flip_time = 5.0

def find_optimal_sln(vessel, target, ref):
	# Get target and vessel current position in sun_rf
	tgt_pos = target.position(ref)
	vessel_pos = vessel.position(ref)
	vessel_vel = vessel.velocity(ref)
	
	# Assume large position error to kick off interration
	position_err = (1e6,1e6,1e6)
	threshold_err = 1e4
	tgt_final_pos = tgt_pos
	now = sc.ut
	
	passnum = 0
	while magnitude(position_err) > threshold_err:
		passnum += 1
		# Get the direction vector pointing towards the target (from vessel to target)
		# dir_vect = tgt_pos - vessel_pos
		dir_vect = vector_subtract(tgt_final_pos, vessel_pos)
		# print('dir_vect = ', dir_vect)
		
		# Get component of vessel velocity in the direction of the target
		v0 = dot_product(dir_vect,vessel_vel)/magnitude(dir_vect)
		vcomp = vector_mult_s(fir_vect, (dot_product(dir_vect,vessel_vel)/dot_product(dir_vect,dir_vect)))
		
		# Get components of vessel velocity perpendicular to the direction of the target
		vperp = vector_subtract(vessel_vel,vcomp)
		
		# normalize the dir_vect
		# dir_vect_n = dir_vect/magnitude(dir_vect)
		
		# Get the distance to the target
		dist = magnitude(dir_vect)
		print('distance = ', dStart)
		
		d3 = (0.5*v0^2)/acc_tgt
		dt3 = v0/acc_tgt
		
		if v0 < 0:
			d2 = (dist-d3)/2
			dt2 = math.sqrt(2*d2/acc_tgt)
			d1 = d2 +  d3
			dt1 = dt2 + dt3
		else:
			d1 = (dist-d3)/2
			dt1 = math.sqrt(2*d1/acc_tgt)
			d2 = d1 +  d3
			dt2 = dt1 + dt3
			
		travelTime = dt1 + dt2 + flip_time
		# General equation for distance traveled under constant acceleration
		# d = v0*t + 0.5a*t^2
		# distance traveled under first part of journey (t0 - t1)
		# t0 = now or start of journey
		# t1 = point of flip
		# t2 = time at end of journey
		# dt1 = t1 - t0 = t1
		# dt2 = t2 - t1
		# v0 = initial velocity
		# v1 = velocity at point of flip
		# v2 = final velocity at end of deceleration burn, generally 0?
		# d1 = v0*dt1 + 0.5a*dt1^2
		# velocity at point of flip (t1)
		# v1 = v0 + a*t1
		# d2 = v1*dt2 - 0.5*a*dt2^2
		# d1 + d2 = dist -> d2 = dist - d1
		# dist - d1 = v1*dt2 - 0.5*a*dt2^2
		# d1 = dist - v1*dt2 + 0.5*a*dt2^2
		# v2 = v1 -a*dt2 = 0
		# v1 = a*dt2
		# Get the time needed to travel that distance via torch ship burns
		# travelTime = 2.0*math.sqrt(dist/acc_tgt)
		
		print('travelTime = ', travelTime)
		
		# Get the target position after the burn
		tgt_new_pos = target.orbit.position_at(now + travelTime, sun_rf)
		# print('tgt_new_pos = ', tgt_new_pos)
		
		# get the position error on this iteration
		position_err = vector_subtract(tgt_new_pos, tgt_final_pos)
		print('%d position_err = %.1f' % (passnum, magnitude(position_err)))
		
		# Update tgt_final_pos
		tgt_final_pos = tgt_new_pos
		# print('tgt_final_pos = ', tgt_final_pos)
		
	sln ={"valid_time": now, "dir_vect": dir_vect, "distance": dStart,
	"travel_time": travelTime, "pos_err": position_err, "passnum": passnum,
	"tgt_pos": tgt_final_pos, "dt1": dt1, "dt2": dt2, "d1": d1, "d2": d2,
	"vperp": vperp, "vcomp": vcomp, "found":True}
	
	# To really make this optimal we should work out when to depart from the current orbit
	# so that our orbital velocity vector at departure is best aligned with dir_vect, then
	# re-optimize at that time just to be sure we're still within tolerance
	
	# Update the trip info
	# update_text1(0.0, travelTime, True)
	# update_text2(dStart, dStart, True)
	
	return sln

# Endless outer loop
while True:
	# Process button press
	if button_clicked():
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
			text_tgt.value = sc.target_body.name
			
			# Require a new optimal solution
			need_optimal_sln = True
		
		# If the active vessel has changed or is not yet defined 
		if vessel != sc.active_vessel:
			# Record the active vessel object

			# Get the vessel auto_pilot for later use
			ap = vessel.auto_pilot
			ap.reference_frame = sun_rf

			# Require a new optimal solution
			need_optimal_sln = True

		# If we need a new optimal solution...
		if need_optimal_sln:
			sln = find_optimal_sln(vessel, target)
			if sln["found"]:
				print('Solution found!')
				print(sln)
				need_optimal_sln = False
				
		# Process button press
		if button_clicked():
			if button.text.content == 'Run':
				if need_optimal_sln
					sln = find_optimal_sln(vessel, target)
					if sln["found"]:
						print('Solution found!')
						print(sln)
						need_optimal_sln = False
				elif: sc.ut - sln["valid_time"] > 100
					need_optimal_sln = True
					sln = find_optimal_sln(vessel, target)
					if sln["found"]:
						print('Solution found!')
						print(sln)
						need_optimal_sln = False
				print("Commencing Run...")
				print("Nulling Verr")
				

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

# set thrust to 0


## Here is where the fun begins!
target = sc.target_body
if target is None:
	print('Assuming target is Mun')
	target = mun
else:
	print('Found target is %s' % target.name)

# Get target and vessel position in sun_rf
tgt_pos = target.position(sun_rf)
vessel_pos = vessel.position(sun_rf)

# Get the direction vector pointing towards the target (from vessel to target)
# dir_vect = tgt_pos - vessel_pos
dir_vect = vector_subtract(tgt_pos, vessel_pos)

# Get the starting distance to the target
dStart = magnitude(dir_vect) - target.sphere_of_influence

travelTime = 2.0*math.sqrt(dStart/(acc_tgt))
startTime = sc.ut
delta_t = sc.ut - startTime

# Update the trip info
update_text1(delta_t, travelTime, True)
update_text2(dStart, dStart, True)

# print('Initial Vessel position (%.1f, %.1f, %.1f)' % vessel_pos)
# print('Initial Mun position (%.1f, %.1f, %.1f)' % tgt_pos)
# print('Initial direction vector (%.1f, %.1f, %.1f)' % dir_vect)
print('Initial distance %.1f km' % (dStart/1000.0))
print('Estimated travel time %.1f seconds' % travelTime)
print('Starting time %.1f' % startTime)
print('Switching SmartASS Off')
# Switch off the MechJeb SmartASS module
sass.autopilot_mode = sass.autopilot_mode.off
sass.update(False)
time.sleep(1)

print('Employing kRPC AP to point to target')
# Make the kRPC autopilot point at the target
ap.target_direction = dir_vect
ap.engage()
time.sleep(1)
if abs(ap.error) > 1:
	ap.wait()
print('Done waiting - should be pointed at target now')
ap.wait()
ap.disengage()

print('Switching SmartASS On in Target+ Mode')
# Configure SmartASS to point towards the target
sass.autopilot_mode = sass.autopilot_mode.target_plus
# Make it so
sass.update(False)
time.sleep(5)

# Accelerate toward target until distance < 1/2 starting distance
# Add while loop here to set and control acceleration
# Hold fixed acceleration towards target until distance is < 1/2
# Cut thrust when while loop ends

# Get the distance to the target
distance = magnitude(dir_vect) - target.sphere_of_influence
delta_t = sc.ut - startTime

print('Starting burn towards target!')
while (distance > dStart/2) and (delta_t < travelTime/2.0):	
	# update the target and vessel position info
	tgt_pos = target.position(sun_rf)
	vessel_pos = vessel.position(sun_rf)
	delta_t = sc.ut - startTime
	
	# Get the direction vector pointing towards the target (from vessel to target)
	# dir_vect = tgt_pos - vessel_pos
	new_dir_vect = vector_subtract(tgt_pos, vessel_pos)

	# Get the remaining distance to the target
	distance = magnitude(new_dir_vect) - target.sphere_of_influence
	delta_t = sc.ut - startTime

	# Set the trust for required acceleration
	acc_max = vessel.available_thrust / vessel.mass / kerbin.surface_gravity
	throttleFraction = acc_tgt_g / acc_max
	# print('Setting thottle to %.3f'% throttleFraction )
	vessel.control.throttle = throttleFraction
	
	# Update the trip info
	update_text1(delta_t, travelTime)
	update_text2(distance, dStart)

# Done with loop
print('Finished Acceleration burn!')
vessel.control.throttle = 0

## Flip vessel
# Get target and vessel position in sun_rf
tgt_pos = target.position(sun_rf)
vessel_pos = vessel.position(sun_rf)

# Get the direction vector pointing away from the target (from target to vessel)
# neg_dir_vect = vessel_pos - tgt_pos
neg_dir_vect = vector_subtract(vessel_pos, tgt_pos)

print('Switching SmartASS Off')
# Switch off the MechJeb SmartASS module
sass.autopilot_mode = sass.autopilot_mode.off
sass.update(False)
time.sleep(1)

print('Employing kRPC AP to point away from target')
# Make the kRPC autopilot point at the target
ap.target_direction = neg_dir_vect
ap.engage()
while abs(ap.error) > 1:
	ap.wait()
print('Done waiting. Should be pointed away from target now')
ap.wait()
ap.disengage()

print('Switching SmartASS On in Target- Mode')
# configure SmartASS to point away from the target
sass.autopilot_mode = sass.autopilot_mode.retrograde
# Make it so
sass.update(False)
time.sleep(1)

# Get target and vessel position in sun_rf
tgt_pos = target.position(sun_rf)
vessel_pos = vessel.position(sun_rf)

# Get the direction vector pointing towards the target (from vessel to target)
# dir_vect = tgt_pos - vessel_pos
dir_vect = vector_subtract(tgt_pos, vessel_pos)

# Get the distance to the target
distance = magnitude(dir_vect)
delta_t = sc.ut - startTime

print('Starting deceleration burn!')
# Decelerate until within target SOI
# Add while loop here to set and control acceleration
# Hold fixed acceleration away from target until within target SOI (or target distance growing!)
# Cut thrust when loop ends
# If within target SOI, then use MJ guidance to get us into an orbit
while distance > 0.75*target.sphere_of_influence:	
	# update the target and vessel position info
	tgt_pos = target.position(sun_rf)
	vessel_pos = vessel.position(sun_rf)
	
	# Get the direction vector pointing towards the target (from vessel to target)
	# dir_vect = vessel_pos - tgt_pos
	new_dir_vect = vector_subtract(vessel_pos, tgt_pos)

	# Get the remaining distance to the target
	distance = magnitude(new_dir_vect)
	delta_t = sc.ut - startTime

	# Set the trust for required acceleration
	acc_max = vessel.available_thrust / vessel.mass / kerbin.surface_gravity
	throttleFraction = acc_tgt_g / acc_max
	# print('Setting throttle to %0.3f' % throttleFraction)
	vessel.control.throttle = throttleFraction
	
	# Update the trip info
	update_text1(delta_t, travelTime)
	update_text2(distance, dStart)

# Done with loop
print('Entered Target SOI!')
vessel.control.throttle = 0

# MJ Set orbital inclination to 0 at fixed time (instantaneous)
# MJ Set periapses to 100 km at fixed time (instantaneous)
# MJ Circularize at periapses
mp = mj.maneuver_planner
setInc = mp.operation_inclination
setPer = mp.operation_periapsis
# setApo = mp.operation_apoapsis
circ = mp.operation_circularize

print('Employing MechJeb to reduce inclination as close to zero as possible right now')
# High thee to a zero-ish inclination orbit!
setInc.new_inclination = 0
setInc.time_selector.lead_time = 0
setInc.time_selector.time_reference = mj.TimeReference.x_from_now
setInc.make_nodes()
execute_nodes()

print('Employing MechJeb to set periapsis to 200km right now')
# High thee to a 100km periapsis orbit!
setPer.new_periapsis = 200000
setPer.time_selector.lead_time = 0
setPer.time_selector.time_reference = mj.TimeReference.x_from_now
setPer.make_nodes()
execute_nodes()

# High thee to a 100km apoapsis orbit!
# setApo.new_periapsis = 100000
# setApo.time_selector.lead_time = 0
# setApo.time_selector.time_reference = mj.TimeReference.x_from_now
# setApo.make_nodes()
# execute_nodes()

print('Employing MechJeb to set circularize the orbit at periapsis')
# Circularize at the next periapsis
circ.time_selector.time_reference = mj.TimeReference.periapsis
circ.make_nodes()
execute_nodes()

print('Employing MechJeb to set circularize the orbit RIGHT NOW DAMN IT!')
# Do it again just in case it wasn't successful the first time
circ.time_selector.lead_time = 0
circ.time_selector.time_reference = mj.TimeReference.x_from_now
circ.make_nodes()
execute_nodes()

print('Done! Enjoy your time the target!')
