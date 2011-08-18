
This code snippet allows you to create an axis aligned
bounding volume tree for a triangle mesh so that you
can do high-speed raycasting.

There are much better implementations of this available
on the internet.  In particular I recommend that you use
OPCODE written by Pierre Terdiman.

@see: http://www.codercorner.com/Opcode.htm

OPCODE does a whole lot more than just raycasting, and
is a rather significant amount of source code.

I am providing this code snippet for the use case where
you *only* want to do quick and dirty optimized raycasting.

I have not done performance testing between this version
and OPCODE; so I don't know how much slower it is.  However,
anytime you switch to using a spatial data structure for
raycasting, you increase your performance by orders and orders
of magnitude; so this implementation should work fine
for simple tools and utilities.

It also serves as a nice sample for people who are trying
to learn the algorithm of how to implement AABB trees.

AABB = Axis Aligned Bounding Volume trees.

http://www.cgal.org/Manual/3.5/doc_html/cgal_manual/AABB_tree/Chapter_main.html


This code snippet was written by John W. Ratcliff on August 18, 2011
and released open source under the MIT. license.

mailto:jratcliffscarab@gmail.com

The official source can be found at:  http://code.google.com/p/raycastmesh/

To run the demo type:

RaycastMesh deer_bound.obj

This will load the Wavefront OBJ file 'deer_bound.obj', then create a 'RaycastMesh'
and then perform one million raycasts against it.

It will then write out an image file 'RaycastMesh.png' to demonstrate that the 
raycasts all worked correctly.  You can feed it other wavefront OBJ files for testing.
