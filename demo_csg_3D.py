from fenics import *
h = 0.25
r = 0.3*h
box = Box(Point(0, 0, 0), Point(1, h, h))
s0 = Sphere(Point(0.3, 0.50*h, 0.50*h), r)
s1 = Sphere(Point(0.5, 0.65*h, 0.65*h), r)
s2 = Sphere(Point(0.7, 0.35*h, 0.35*h), r)
domain = box - s0 - s1 - s2
 
# Generate mesh
mesh = generate_mesh(domain, 32)
