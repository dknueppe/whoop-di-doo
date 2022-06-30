#!/usr/bin/env python

# tiny whoop frame. design currently makes the impression of being configurable
# but is barely... changing certain values will likely break the design as of now.

import cadquery as cq
from math import sin, cos, tau
import numpy as np
import numpy.linalg as la

# next part stolen from here https://github.com/CadQuery/cadquery/issues/371

from cadquery import Workplane
from typing import Dict
from math import pi
from cadquery.occ_impl.shapes import TopAbs_Orientation, Shell, Edge


def edge_angle_map(shell: Shell, types=["CIRCLE", "LINE"]) -> Dict[Edge, float]:
    """returns a dictionary where the keys are edges and the values are angles
    between the adjoining faces, with negative interior angles and positive
    exterior angles

    Note that angles are not generally well defined for edges other than
    circles and lines. It may be well defined for some instances of other
    edge types depending on their construction.  This could be tested for
    heuristically, but for now I'm only returning edges for lines and
    circles by default.
    """
    if not shell.Closed():
        raise RuntimeError("Shell should be closed")
    d = shell._entitiesFrom("Edge", "Face")
    # seams in sphere's and cylinders only touch one face.  Also see note above:
    d = dict((k, v) for k, v in d.items() if len(v) == 2 and k.geomType() in types)
    out = {}
    for e, (f0, f1) in d.items():
        pt = e.positionAt(0)
        v0 = f0.normalAt(pt)
        v1 = f1.normalAt(pt)
        a = 180 * v0.getAngle(v1) / pi
        n = e.tangentAt(0)
        det = (
            n.x * (v0.y * v1.z - v0.z * v1.y)
            - n.y * (v0.x * v1.z - v0.z * v1.x)
            + n.z * (v0.x * v1.y - v0.y * v1.x)
        )
        if e.wrapped.Orientation() != TopAbs_Orientation.TopAbs_FORWARD:
            det *= -1
        out[e] = -a if det < 0 else a
    return out

def inside_edges(x: Workplane) -> list[Edge]:
    """select the edges with negative angles between the faces"""
    mappings = [edge_angle_map(s) for s in x.shells().objects if s.Closed()]
    edges = [[k for k, v in d.items() if v < 0] for d in mappings]
    return [e for el in edges for e in el]

def outside_edges(x: Workplane) -> list[Edge]:
    """select the edges with negative angles between the faces"""
    mappings = [edge_angle_map(s) for s in x.shells().objects if s.Closed()]
    edges = [[k for k, v in d.items() if v < 0] for d in mappings]
    return [e for el in edges for e in el]

def select_edges_by_angle(x: Workplane, min=-180, max=180) -> list[Edge]:
    """select the edges with negative angles between the faces"""
    mappings = [edge_angle_map(s) for s in x.shells().objects if s.Closed()]
    edges = [[k for k, v in d.items() if min < v < max] for d in mappings]
    return [e for el in edges for e in el]

############ up until here https://github.com/CadQuery/cadquery/issues/371

mm_per_inch = 25.4
mat_thicc = 2
fillet_rad = 2

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
        .circle(4.1)
        .extrude(mat_thicc) )

props = ( cq.Workplane("front")
        .pushPoints(motor_pos)
        .circle(prop_size / 2)
        .extrude(4)
        .translate((0,0,15)))

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
                .trapezoid(2.85, float(length), 88.5)
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

# Loads of whilly nilly magic numbers and fiddling ahead, lost all
# interest in making things 'parametric' and just wanted this one piece done.
# Should refactor later.
cam_mnt = ( cq.Workplane("front")
        .sketch()
        .trapezoid(20, 10, 75)
        .finalize()
        .extrude(mat_thicc)
        .translate((0, fc_center_dist + 8))
        )
cam_mnt = cam_mnt.union( cq.Workplane("front")
        .sketch()
        .trapezoid(20, -3, 150)
        .finalize()
        .extrude(mat_thicc)
        .translate((0, fc_center_dist + 2.5))
        )
c_mnt_cutout = ( cq.Workplane("front")
        .rect(15, 2.5) 
        .extrude(mat_thicc) )

c_mnt_cutout = c_mnt_cutout.union( cq.Workplane("front")
        .rect(15 - cos(1/8 * tau), 2.5 - cos(1/8 * tau), forConstruction=True)
        .vertices()
        .circle(0.5)
        .extrude(mat_thicc))

bow_tie = (cq.Workplane("front")
        .sketch()
        .trapezoid(3.3, 8, 80)
        .finalize()
        .extrude(mat_thicc)
        .rotate((0,0,0), (0,0,1), 90)
        .translate((4,0)))

bow_tie = bow_tie.union(bow_tie.mirror(mirrorPlane="ZY", basePointVector=(0, 0 ,0)))

frame = frame.union(cam_mnt)
frame = frame.union(cq.Workplane("front")
        .rect(1.5, 2 * fc_center_dist)
        .extrude(mat_thicc))

frame = frame.union(cq.Workplane("front")
        .rect(7, 10)
        .extrude(mat_thicc)
        .translate((0, -fc_center_dist -10/2)))

frame = frame.union(bow_tie.translate((0, -5)))
frame = frame.union(bow_tie.translate((0, 5)))


zip_cutout = ( cq.Workplane("front")
        .rect(4.7, 2) 
        .extrude(mat_thicc) )

zip_cutout = zip_cutout.union( cq.Workplane("front")
        .rect(4.7 - cos(1/8 * tau), 2 - cos(1/8 * tau), forConstruction=True)
        .vertices()
        .circle(0.5)
        .extrude(mat_thicc))

frame = frame.cut(c_mnt_cutout.translate((7.5, fc_center_dist + 10)))
frame = frame.cut(zip_cutout.translate((0, -fc_center_dist - 7.5)))

frame = frame.cut(cq.Workplane("front")
        .ellipse(7,2.5)
        .extrude(mat_thicc)
        .translate((0,fc_center_dist + 5)))

for i , motor in enumerate(motor_pos):
    frame = frame.union( cq.Workplane("front")
            .polygon(3, 14)
            .extrude(mat_thicc)
            .edges("|Z").fillet(2)
            .rotateAboutCenter((0,0,1), i * 90 +180 -135)
            .translate((*motor, 0)))

frame = frame.newObject(inside_edges(frame)).fillet(fillet_rad)

for i , motor in enumerate(motor_pos):
    frame = (frame.cut(fc_mount)
                .cut(motor_mount.rotateAboutCenter((0,0,1), i * 90 -45)
                .translate((*motor, 0))))

### Camera mount

front_plate = (cq.Workplane("front")
        .moveTo(0, 25)
        .rect(10, 50)
        .moveTo(0, -3.75)
        .rect(15.5, 7.5)
        .moveTo(0, 23)
        .circle(6.5)
        .moveTo(4, 47)
        .rect(26, 15)
        .moveTo(3.5, 56.55)
        .rect(15, 4.1)
        .extrude(2)
        .moveTo(0, 23)
        .circle(4.95)
        .moveTo(7, 47)
        .rect(17, 12)
        .moveTo(-9, 47)
        .circle(5)
        .moveTo(0, 19)
        .cutBlind(2)
        .moveTo(7, 59.6)
        .rect(50, 2)
        .extrude(2)
        .moveTo(0, 18)
        .rect(6, 34)
        .moveTo(6.5, 56.5)
        .circle(2.5)
        .cutBlind(2))

front_plate = front_plate.newObject(inside_edges(front_plate)).fillet(fillet_rad)
front_plate = front_plate.cut(zip_cutout.translate((0,56.5)))
front_plate = front_plate.cut(zip_cutout.translate((0,38)))
front_plate = front_plate.cut(c_mnt_cutout.translate((7.5, -3.5)))

cq.exporters.export(frame, "base_plate.step")
frame = frame.add(front_plate.rotate(origin, (0,1,0), 180).rotate(origin, (1,0,0), 105).translate((0, 26, 4.2)))
cq.exporters.export(frame, "frame.step")
cq.exporters.export(front_plate, "front_plate.step")

