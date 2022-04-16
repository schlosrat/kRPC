# kRPC
KRPC python scripts

This is a collection of Python3 scripts for use with kRPC and kRPC.MechJeb. These scripts have been developed using the following SW versions.

Python 3.6.8
kRPC 0.4.9 compiled for KSP 1.12.1
kRPC.MechJeb 0.6.1 (built for MechJeb 2.8.3.0)
KSP 1.12.3
MechJeb 2.14.1

NOTE: Due to changes in MechJeb between 2.8 and 2.14 there are some MechJeb features that are not accessible via kRPC.MechJeb 0.6.1. nevertheless, 0.6.1 is the most recent version and is able to control a number of things in MechJeb 2.14 which has been sufficent for my purposes.

This project began as a means to meet my need for a Brachristochrone (minimum time) trajectory control process useful with torch ships able to accelerate continuously. Numerous KSP mods provide suitable engines with extremely high ISP and morerate to high thrust. Such a trajectory is not what you would get with a low thrust / high ISP engine like a solar sail or ion engine. To accomplish this sort of trajectory the vessel must be capable of prolonged thrust at an thrust greater than the force of the Sun's gravity.
