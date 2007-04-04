#!/usr/bin/env python
import sys
from math import *
from euclid import *



def equatorialToGalactic():
    r1 = Matrix4.new_rotatez(-282.8592*pi/180)
    r2 = Matrix4.new_rotatex(-62.8717*pi/180)
    r3 = Matrix4.new_rotatez(32.93224*pi/180)
    
    return r3*r2*r1

def galacticToEquatorial():
    r1 = Matrix4.new_rotatez(+282.8592*pi/180)
    r2 = Matrix4.new_rotatex(+62.8717*pi/180)
    r3 = Matrix4.new_rotatez(-32.93224*pi/180)
    
    return r1*r2*r3

def Cel2Gal(ra,dec):
    ra_rad  = float(ra)*pi/180
    dec_rad = float(dec)*pi/180

    v1 = Vector3(cos(ra_rad)*cos(dec_rad ), sin(ra_rad)*cos(dec_rad ), sin(dec_rad ))
    v2 = equatorialToGalactic()*v1
    l = atan2(v2[1],v2[0])*180/pi
    if l < 0:
	l=l+360
    b = asin(v2[2])*180/pi
    print 'ra= ',ra,' dec= ',dec
    print 'l= ',l,' b= ',b
    return [l,b]

def Gal2Cel(l,b):
    
    l_rad = float(l)*pi/180
    b_rad = float(b)*pi/180
    v1 = Vector3(cos(l_rad)*cos(b_rad ), sin(l_rad)*cos(b_rad ), sin(b_rad ))
    v2 = galacticToEquatorial()*v1
    ra = atan2(v2[1],v2[0])*180/pi
    if ra < 0:
	ra=ra+360
    dec = asin(v2[2])*180/pi
    print 'ra(deg)= ',ra,' ra(hr) = ',ra*24./360,', dec= ',dec
    print 'l= ',l,' b= ',b
    return [ra,dec]


if __name__ == '__main__':
    filename=sys.argv[1]
    
    fin  = file(filename,'r')
    fout = file('l_b_coords.txt','w')
    
    L=fin.readlines()
    for l in L:
	line = l.split()
	time=line[0]
	ra=line[1]
	dec=line[2]
	l_b=Cel2Gal(ra,dec)
	l=l_b[0]
	b=l_b[1]
	newline = str(time)+' %.4f %.4f ' % (l , b) + ' \n '
	fout.writelines(newline)
