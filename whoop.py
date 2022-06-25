#!/usr/bin/env python

# tiny whoop frame. design currently makes the impression of being configurable
# but is barely... changing certain values will likely break the design as of now.

import cadquery as cq
from math import sin, cos, tau
import numpy as np
import numpy.linalg as la

mm_per_inch = 25.4
mat_thicc = 2
fillet_rad = 3

# generate "n_holes" amount of points, that are equally spaced on a circle of
# radius "mot_mnt_dia", placing the first one on the Y-axis facing "up".
n_holes = 3
mot_mnt_dia = 6.6/2
holes = np.array([ ( mot_mnt_dia * cos((n / n_holes + 1/4) * tau), 
                     mot_mnt_dia * sin((n / n_holes + 1/4) * tau) )
                     for n in range(n_holes) ])
motor_mount = ( cq.Workplane("front")
        .circle(2)
        .pushPoints(holes)
        .circle(1.4 / 2)
        .extrude(mat_thicc) )

fc_mount_dist = 26 # non square flight controller mount patterns are not supported
fc_center_dist = (2 * fc_mount_dist ** 2) ** 0.5 / 2
# holes are arranged in a diamond pattern and ordered counter clockwise
# The first hole corresponds to the "3" position on the clock
fc_holes = np.array([ ( fc_center_dist * cos((n / 4) * tau), 
                        fc_center_dist * sin((n / 4) * tau) )
                        for n in range(4) ])
fc_mount = ( cq.Workplane("front")
        .pushPoints(fc_holes)
        .circle(0.5)
        .extrude(mat_thicc) )

prop_size = 1.6 * mm_per_inch
safety_dist = prop_size * 0.3
prop_dist = round(prop_size + safety_dist)
stretch = 1
deadcat = 0
wheelbase_rad = ( 2 * (prop_dist / 2) **2 ) ** 0.5
n_motors = 4
motor_pos = np.array([ ( wheelbase_rad * cos((n / n_motors + 1/8) * tau),
                         wheelbase_rad * sin((n / n_motors + 1/8) * tau) )
                         for n in range(n_motors) ])
motor_plates = ( cq.Workplane("front")
        .pushPoints(motor_pos)
        .circle(4.25)
        .extrude(mat_thicc) )

frame = motor_plates
#calculate lengths and angle of trusses for arms
origin = (0,0,0)
unity_ref = np.array([ ( cos((n / n_motors) * tau),
                         sin((n / n_motors) * tau) )
                         for n in range(n_motors) ])

for j in range(2):
    for k, (ref, motor) in enumerate(zip(unity_ref, motor_pos)):
        i = (k+j) % len(fc_holes)
        truss_vec = motor - fc_holes[i]
        length = la.norm(truss_vec)
        inner = np.inner(truss_vec, ref)
        norms = la.norm(truss_vec) * la.norm(ref)
        angle = np.arccos(np.clip(inner / norms, -1.0, 1.0))
        deg = np.rad2deg(angle) + (i - n_motors / 4 - j) * 360 / n_motors
        
        truss = ( cq.Workplane("front")
                .sketch()
                .trapezoid(2.85, float(length), 88)
                .finalize()
                .extrude(mat_thicc)
                .translate((fc_holes[i][0],float(length) / 2 + fc_holes[i][1]))
                .rotate((*fc_holes[i], 0), (*fc_holes[i], 1), deg))
        frame = frame.union(truss)

# inner support for flight controller
frame = frame.union( cq.Workplane("front")
        .rect(fc_mount_dist + 1, fc_mount_dist + 1)
        .extrude(mat_thicc)
        .rotateAboutCenter((0,0,1), 45) )
frame = frame.cut( cq.Workplane("front")
        .rect(fc_mount_dist - 1, fc_mount_dist -1)
        .extrude(mat_thicc)
        .rotateAboutCenter((0,0,1), 45) )

frame = frame.edges("|Z").fillet(2)

#frame = frame.edges(">Z").fillet(0.2)
#frame = frame.edges("<Z").fillet(0.2)

for i , motor in enumerate(motor_pos):
    frame = (frame.cut(fc_mount)
                .cut(motor_mount.rotateAboutCenter((0,0,1), i * 90 -45)
                .translate((*motor, 0))))
    
cq.exporters.export(frame, "frame.step")

