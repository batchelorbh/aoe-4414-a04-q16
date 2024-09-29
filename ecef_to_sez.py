#!/usr/bin/env python
# ecef_to_sez.py
#
# Converts ECEF reference frame components to ECEF
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
#
# Written by Blake Batchelor, batchelorbh@vt.edu
# Other contributors: None
#
# Parameters:
#    o_x_km              X component of EFEF origin of SEZ frame in km
#    o_y_km              Y component of EFEF origin of SEZ frame in km
#    o_z_km              Z component of EFEF origin of SEZ frame in km
#    x_km                X component of the ECEF position in km
#    y_km                Y component of the ECEF position in km
#    z_km                Z component of the ECEF position in km
#
# Output:
#    Prints SEZ vector components in km
#
# Revision history:
#    09/28/2024          Script created
#
###############################################################################

#Import relevant modules
import sys
from math import pi, sqrt, sin, cos, atan2, isnan, asin, atan

#Define constants
R_E_KM = 6378.137 #km
E_E = 0.081819221456
RAD_TO_DEG = 180.0 / pi

#Calculate denominator for SE and CE equations
def calc_denom(ecc, lat_rad):
   return sqrt(1 - ecc**2 * sin(lat_rad)**2)

#Pre-initialize input parameters
o_x_km = float('nan') #x component of EFEF origin of SEZ frame in km
o_y_km = float('nan') #y component of EFEF origin of SEZ frame in km
o_z_km = float('nan') #z component of EFEF origin of SEZ frame in km
x_km = float('nan') #x component of the ECEF position in km
y_km = float('nan') #y component of the ECEF position in km
z_km = float('nan') #z component of the ECEF position in km

#Arguments are strings by default
if len(sys.argv) == 7:
   o_x_km = float(sys.argv[1])
   o_y_km = float(sys.argv[2])
   o_z_km = float(sys.argv[3])
   x_km = float(sys.argv[4])
   y_km = float(sys.argv[5])
   z_km = float(sys.argv[6])
else:
   print('Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km')
   sys.exit()

#Main body of script

#Get vector components
x_vec_km = x_km - o_x_km
y_vec_km = y_km - o_y_km
z_vec_km = z_km - o_z_km

#Calculate longitude
lon_rad = atan2(o_y_km, o_x_km)
lon_deg = lon_rad * RAD_TO_DEG

#Initialize lat_rad, r_lon_km, r_z_km
lat_rad = asin(z_vec_km / sqrt(o_x_km**2 + o_y_km**2 + o_z_km**2))
r_lon_km = sqrt(o_x_km**2 + o_y_km**2)
prev_lat_rad = float('nan')

#Iteratively find latitude
c_E = float('nan')
count = 0
while (isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad) > 10e-7) and count < 5:
  denom = calc_denom(E_E, lat_rad)
  c_E = R_E_KM / denom
  prev_lat_rad = lat_rad
  lat_rad = atan((o_z_km + c_E * (E_E**2) * sin(lat_rad)) / r_lon_km)
  count = count + 1

#Calculate HAE
hae_km = r_lon_km / cos(lat_rad) - c_E

#Rotation matrices
s_km = x_vec_km * cos(lon_rad) + y_vec_km * sin(lon_rad) \
       * sin(lat_rad) - z_vec_km * cos(lat_rad)

e_km = x_vec_km * -sin(lon_rad) + y_vec_km * cos(lon_rad)

z_km = (x_vec_km * cos(lon_rad) + y_vec_km * sin(lon_rad)) \
       * cos(lat_rad) + z_vec_km * sin(lat_rad)

#Print SEZ components
print(s_km)
print(e_km)
print(z_km)

'''
print(lat_rad * RAD_TO_DEG)
print(lon_deg)
print(hae_km)
'''
