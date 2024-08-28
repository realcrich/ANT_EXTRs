# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 11:07:01 2024

@author: Collin
"""

# written by Tyler Sutterley, just included in this library for ease of use for user

def map_ll(alon, alat, HEM='S', FORMAT='dict'):
	import numpy as np

	rad_e = 6378137.0#-- WGS84
	ecc2 = 0.00669437999015#-- WGS84
	ecc = np.sqrt(ecc2)

	#-- Standard parallel - latitude with no distortion = +70/-71.
	#-- slat: unsigned latitude of origin (N: 70, S: 71)
	#-- sn: sign of latitude (N: +1, S: -1)
	#-- xlam: longitude of center (N: 45, S: 0)
	if HEM in ('N','n'):#-- Northern Hemisphere
		slat = 70.0
		sn = 1.0
		xlam = 45.0
	elif HEM in ('S','s'):#-- Southern Hemisphere
		slat = 71.0
		sn = -1.0
		xlam = 0.0
	elif (HEM == 'GDEM'):#-- Bamber Greenland 5km DEM
		slat = 71.0
		sn = 1.0
		xlam = 39.0
	else:
		slat,sn,xlam = HEM

	alat = sn*alat
	alon = sn*alon
	rlat = alat

	t1 = np.tan(np.pi/4.0 - rlat*np.pi/(2.0*180.0)) / \
		((1.0 - ecc*np.sin(rlat*np.pi/180.0))/ \
		(1.0 + ecc*np.sin(rlat*np.pi/180.0)))**(ecc/2.0)

	t2 = np.tan(np.pi/4.0 - slat*np.pi/(2.0*180.0)) / \
		((1.0 - ecc*np.sin(slat*np.pi/180.0)) / \
		(1.0 + ecc*np.sin(slat*np.pi/180.0)))**(ecc/2.0)
	#-- m at standard latitude
	cm = np.cos(slat*np.pi/180.0) / \
		np.sqrt(1.0-ecc2*(np.sin(slat*np.pi/180.0)**2.0))
	#-- radius of latitude circle
	rho = rad_e*cm*t1/t2
	#-- polar stereographic x and y
	x = ( rho*sn*np.sin((alon+xlam)*np.pi/180.0))
	y = (-rho*sn*np.cos((alon+xlam)*np.pi/180.0))

	#-- return coordinates in output format (default python dictionary)
	if (FORMAT == 'dict'):
		return dict(x=x,y=y)
	elif (FORMAT == 'tuple'):
		return (x,y)
	elif (FORMAT == 'zip'):
		return zip(x,y)