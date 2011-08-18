#ifndef RAYCAST_MESH_H

#define RAYCAST_MESH_H

typedef float RmReal;
typedef unsigned int RmUint32;

class RaycastMesh
{
public:
	virtual bool raycast(const RmReal *from,const RmReal *to,RmReal *hitLocation,RmReal *hitNormal,RmReal *hitDistance) = 0;
	virtual const RmReal * getBoundMin(void) const = 0; // return the minimum bounding box
	virtual const RmReal * getBoundMax(void) const = 0; // return the maximum bounding box.
	virtual void release(void) = 0;
protected:
	virtual ~RaycastMesh(void) { };
};


RaycastMesh * createRaycastMesh(RmUint32 vcount,		// The number of vertices in the source triangle mesh
								const RmReal *vertices,		// The array of vertex positions in the format x1,y1,z1..x2,y2,z2.. etc.
								RmUint32 tcount,		// The number of triangles in the source triangle mesh
								const RmUint32 *indices, // The triangle indices in the format of i1,i2,i3 ... i4,i5,i6, ...
								RmUint32 maxDepth=15,	// Maximum recursion depth for the triangle mesh.
								RmUint32 minLeafSize=4,	// minimum triangles to treat as a 'leaf' node.
								RmReal	minAxisSize=0.01f	// once a particular axis is less than this size, stop sub-dividing.
								);


#endif
