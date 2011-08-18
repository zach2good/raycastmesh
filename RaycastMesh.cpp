#include "RaycastMesh.h"
#include <math.h>
#include <assert.h>
#include <vector>

// This code snippet allows you to create an axis aligned bounding volume tree for a triangle mesh so that you can do
// high-speed raycasting.
//
// There are much better implementations of this available on the internet.  In particular I recommend that you use 
// OPCODE written by Pierre Terdiman.
// @see: http://www.codercorner.com/Opcode.htm
//
// OPCODE does a whole lot more than just raycasting, and is a rather significant amount of source code.
//
// I am providing this code snippet for the use case where you *only* want to do quick and dirty optimized raycasting.
// I have not done performance testing between this version and OPCODE; so I don't know how much slower it is.  However,
// anytime you switch to using a spatial data structure for raycasting, you increase your performance by orders and orders 
// of magnitude; so this implementation should work fine for simple tools and utilities.
//
// It also serves as a nice sample for people who are trying to learn the algorithm of how to implement AABB trees.
// AABB = Axis Aligned Bounding Volume trees.
//
// http://www.cgal.org/Manual/3.5/doc_html/cgal_manual/AABB_tree/Chapter_main.html
//
//
// This code snippet was written by John W. Ratcliff on August 18, 2011 and released under the MIT. license.
//
// mailto:jratcliffscarab@gmail.com
//
// The official source can be found at:  http://code.google.com/p/raycastmesh/
//
// 

#pragma warning(disable:4100)

#define USE_BRUTE_FORCE 0

namespace RAYCAST_MESH
{

typedef std::vector< RmUint32 > TriVector;

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/**
		*	A method to compute a ray-AABB intersection.
		*	Original code by Andrew Woo, from "Graphics Gems", Academic Press, 1990
		*	Optimized code by Pierre Terdiman, 2000 (~20-30% faster on my Celeron 500)
		*	Epsilon value added by Klaus Hartmann. (discarding it saves a few cycles only)
		*
		*	Hence this version is faster as well as more robust than the original one.
		*
		*	Should work provided:
		*	1) the integer representation of 0.0f is 0x00000000
		*	2) the sign bit of the RmReal is the most significant one
		*
		*	Report bugs: p.terdiman@codercorner.com
		*
		*	\param		aabb		[in] the axis-aligned bounding box
		*	\param		origin		[in] ray origin
		*	\param		dir			[in] ray direction
		*	\param		coord		[out] impact coordinates
		*	\return		true if ray intersects AABB
		*/
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		#define RAYAABB_EPSILON 0.00001f
		//! Integer representation of a RmRealing-point value.
		#define IR(x)	((RmUint32&)x)

		bool intersectRayAABB(const RmReal MinB[3],const RmReal MaxB[3],const RmReal origin[3],const RmReal dir[3],RmReal coord[3])
		{
			bool Inside = true;
			RmReal MaxT[3];
			MaxT[0]=MaxT[1]=MaxT[2]=-1.0f;

			// Find candidate planes.
			for(RmUint32 i=0;i<3;i++)
			{
				if(origin[i] < MinB[i])
				{
					coord[i]	= MinB[i];
					Inside		= false;

					// Calculate T distances to candidate planes
					if(IR(dir[i]))	MaxT[i] = (MinB[i] - origin[i]) / dir[i];
				}
				else if(origin[i] > MaxB[i])
				{
					coord[i]	= MaxB[i];
					Inside		= false;

					// Calculate T distances to candidate planes
					if(IR(dir[i]))	MaxT[i] = (MaxB[i] - origin[i]) / dir[i];
				}
			}

			// Ray origin inside bounding box
			if(Inside)
			{
				coord[0] = origin[0];
				coord[1] = origin[1];
				coord[2] = origin[2];
				return true;
			}

			// Get largest of the maxT's for final choice of intersection
			RmUint32 WhichPlane = 0;
			if(MaxT[1] > MaxT[WhichPlane])	WhichPlane = 1;
			if(MaxT[2] > MaxT[WhichPlane])	WhichPlane = 2;

			// Check final candidate actually inside box
			if(IR(MaxT[WhichPlane])&0x80000000) return false;

			for(RmUint32 i=0;i<3;i++)
			{
				if(i!=WhichPlane)
				{
					coord[i] = origin[i] + MaxT[WhichPlane] * dir[i];
					#ifdef RAYAABB_EPSILON
					if(coord[i] < MinB[i] - RAYAABB_EPSILON || coord[i] > MaxB[i] + RAYAABB_EPSILON)	return false;
					#else
					if(coord[i] < MinB[i] || coord[i] > MaxB[i])	return false;
					#endif
				}
			}
			return true;	// ray hits box
		}




		bool intersectLineSegmentAABB(const RmReal bmin[3],const RmReal bmax[3],const RmReal p1[3],const RmReal dir[3],RmReal &dist,RmReal intersect[3])
		{
			bool ret = false;

			if ( dist > RAYAABB_EPSILON )
			{
				ret = intersectRayAABB(bmin,bmax,p1,dir,intersect);
				if ( ret )
				{
					RmReal dx = p1[0]-intersect[0];
					RmReal dy = p1[1]-intersect[1];
					RmReal dz = p1[2]-intersect[2];
					RmReal d = dx*dx+dy*dy+dz*dz;
					if ( d < dist*dist )
					{
						dist = sqrtf(d);
					}
					else
					{
						ret = false;
					}
				}
			}
			return ret;
		}




/* a = b - c */
#define vector(a,b,c) \
	(a)[0] = (b)[0] - (c)[0];	\
	(a)[1] = (b)[1] - (c)[1];	\
	(a)[2] = (b)[2] - (c)[2];



#define innerProduct(v,q) \
	((v)[0] * (q)[0] + \
	(v)[1] * (q)[1] + \
	(v)[2] * (q)[2])

#define crossProduct(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];


static inline bool rayIntersectsTriangle(const RmReal *p,const RmReal *d,const RmReal *v0,const RmReal *v1,const RmReal *v2,RmReal &t)
{
	RmReal e1[3],e2[3],h[3],s[3],q[3];
	RmReal a,f,u,v;

	vector(e1,v1,v0);
	vector(e2,v2,v0);
	crossProduct(h,d,e2);
	a = innerProduct(e1,h);

	if (a > -0.00001 && a < 0.00001)
		return(false);

	f = 1/a;
	vector(s,p,v0);
	u = f * (innerProduct(s,h));

	if (u < 0.0 || u > 1.0)
		return(false);

	crossProduct(q,s,e1);
	v = f * innerProduct(d,q);
	if (v < 0.0 || u + v > 1.0)
		return(false);
	// at this stage we can compute t to find out where
	// the intersection point is on the line
	t = f * innerProduct(e2,q);
	if (t > 0) // ray intersection
		return(true);
	else // this means that there is a line intersection
		// but not a ray intersection
		return (false);
}

static RmReal computePlane(const RmReal *A,const RmReal *B,const RmReal *C,RmReal *n) // returns D
{
	RmReal vx = (B[0] - C[0]);
	RmReal vy = (B[1] - C[1]);
	RmReal vz = (B[2] - C[2]);

	RmReal wx = (A[0] - B[0]);
	RmReal wy = (A[1] - B[1]);
	RmReal wz = (A[2] - B[2]);

	RmReal vw_x = vy * wz - vz * wy;
	RmReal vw_y = vz * wx - vx * wz;
	RmReal vw_z = vx * wy - vy * wx;

	RmReal mag = sqrt((vw_x * vw_x) + (vw_y * vw_y) + (vw_z * vw_z));

	if ( mag < 0.000001f )
	{
		mag = 0;
	}
	else
	{
		mag = 1.0f/mag;
	}

	RmReal x = vw_x * mag;
	RmReal y = vw_y * mag;
	RmReal z = vw_z * mag;


	RmReal D = 0.0f - ((x*A[0])+(y*A[1])+(z*A[2]));

	n[0] = x;
	n[1] = y;
	n[2] = z;

	return D;
}


/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-Möller                              */
/* Function: int triBoxOverlap(RmReal boxcenter[3],      */
/*          RmReal boxhalfsize[3],RmReal triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/
#include <math.h>
#include <stdio.h>

#define X 0
#define Y 1
#define Z 2

#define CROSS(dest,v1,v2) \
	dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
	dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
	dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2) \
	dest[0]=v1[0]-v2[0]; \
	dest[1]=v1[1]-v2[1]; \
	dest[2]=v1[2]-v2[2];

#define FINDMINMAX(x0,x1,x2,min,max) \
	min = max = x0;   \
	if(x1<min) min=x1;\
	if(x1>max) max=x1;\
	if(x2<min) min=x2;\
	if(x2>max) max=x2;

int planeBoxOverlap(RmReal normal[3],RmReal d, RmReal maxbox[3])
{
	int q;
	RmReal vmin[3],vmax[3];
	for(q=X;q<=Z;q++)
	{
		if(normal[q]>0.0f)
		{
			vmin[q]=-maxbox[q];
			vmax[q]=maxbox[q];
		}
		else
		{
			vmin[q]=maxbox[q];
			vmax[q]=-maxbox[q];
		}
	}
	if(DOT(normal,vmin)+d>0.0f) return 0;
	if(DOT(normal,vmax)+d>=0.0f) return 1;

	return 0;
}


/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)             \
	p0 = a*v0[Y] - b*v0[Z];                    \
	p2 = a*v2[Y] - b*v2[Z];                    \
	if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;

#define AXISTEST_X2(a, b, fa, fb)              \
	p0 = a*v0[Y] - b*v0[Z];                    \
	p1 = a*v1[Y] - b*v1[Z];                    \
	if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)             \
	p0 = -a*v0[X] + b*v0[Z];                   \
	p2 = -a*v2[X] + b*v2[Z];                       \
	if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;

#define AXISTEST_Y1(a, b, fa, fb)              \
	p0 = -a*v0[X] + b*v0[Z];                   \
	p1 = -a*v1[X] + b*v1[Z];                       \
	if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)             \
	p1 = a*v1[X] - b*v1[Y];                    \
	p2 = a*v2[X] - b*v2[Y];                    \
	if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return 0;

#define AXISTEST_Z0(a, b, fa, fb)              \
	p0 = a*v0[X] - b*v0[Y];                \
	p1 = a*v1[X] - b*v1[Y];                    \
	if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return 0;

int triBoxOverlap(RmReal boxcenter[3],RmReal boxhalfsize[3],RmReal triverts[3][3])
{

	/*    use separating axis theorem to test overlap between triangle and box */
	/*    need to test for overlap in these directions: */
	/*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
	/*       we do not even need to test these) */
	/*    2) normal of the triangle */
	/*    3) crossproduct(edge from tri, {x,y,z}-directin) */
	/*       this gives 3x3=9 more tests */
	RmReal v0[3],v1[3],v2[3];
	RmReal min,max,d,p0,p1,p2,rad,fex,fey,fez;
	RmReal normal[3],e0[3],e1[3],e2[3];

	/* This is the fastest branch on Sun */
	/* move everything so that the boxcenter is in (0,0,0) */
	SUB(v0,triverts[0],boxcenter);
	SUB(v1,triverts[1],boxcenter);
	SUB(v2,triverts[2],boxcenter);

	/* compute triangle edges */
	SUB(e0,v1,v0);      /* tri edge 0 */
	SUB(e1,v2,v1);      /* tri edge 1 */
	SUB(e2,v0,v2);      /* tri edge 2 */

	/* Bullet 3:  */
	/*  test the 9 tests first (this was faster) */
	fex = fabs(e0[X]);
	fey = fabs(e0[Y]);
	fez = fabs(e0[Z]);
	AXISTEST_X01(e0[Z], e0[Y], fez, fey);
	AXISTEST_Y02(e0[Z], e0[X], fez, fex);
	AXISTEST_Z12(e0[Y], e0[X], fey, fex);

	fex = fabs(e1[X]);
	fey = fabs(e1[Y]);
	fez = fabs(e1[Z]);
	AXISTEST_X01(e1[Z], e1[Y], fez, fey);
	AXISTEST_Y02(e1[Z], e1[X], fez, fex);
	AXISTEST_Z0(e1[Y], e1[X], fey, fex);

	fex = fabs(e2[X]);
	fey = fabs(e2[Y]);
	fez = fabs(e2[Z]);
	AXISTEST_X2(e2[Z], e2[Y], fez, fey);
	AXISTEST_Y1(e2[Z], e2[X], fez, fex);
	AXISTEST_Z12(e2[Y], e2[X], fey, fex);

	/* Bullet 1: */
	/*  first test overlap in the {x,y,z}-directions */
	/*  find min, max of the triangle each direction, and test for overlap in */
	/*  that direction -- this is equivalent to testing a minimal AABB around */
	/*  the triangle against the AABB */

	/* test in X-direction */
	FINDMINMAX(v0[X],v1[X],v2[X],min,max);
	if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return 0;

	/* test in Y-direction */
	FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
	if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return 0;

	/* test in Z-direction */
	FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
	if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return 0;

	/* Bullet 2: */
	/*  test if the box intersects the plane of the triangle */
	/*  compute plane equation of triangle: normal*x+d=0 */
	CROSS(normal,e0,e1);
	d=-DOT(normal,v0);  /* plane eq: normal.x+d=0 */
	if(!planeBoxOverlap(normal,d,boxhalfsize)) return 0;

	return 1;   /* box and triangle overlaps */
}

#define TRI_EOF 0xFFFFFFFF

enum AxisAABB
{
	AABB_XAXIS,
	AABB_YAXIS,
	AABB_ZAXIS
};

enum ClipCode
{
	OLEFT   =       (1<<0),
	ORIGHT  =       (1<<1),
	OTOP    =       (1<<2),
	OBOTTOM =       (1<<3),
	OFRONT  =       (1<<4),
	OBACK   =       (1<<5),
};


class BoundsAABB
{
public:


	void setMin(const RmReal *v)
	{
		mMin[0] = v[0];
		mMin[1] = v[1];
		mMin[2] = v[2];
	}

	void setMax(const RmReal *v)
	{
		mMax[0] = v[0];
		mMax[1] = v[1];
		mMax[2] = v[2];
	}

	void setMin(RmReal x,RmReal y,RmReal z)
	{
		mMin[0] = x;
		mMin[1] = y;
		mMin[2] = z;
	}

	void setMax(RmReal x,RmReal y,RmReal z)
	{
		mMax[0] = x;
		mMax[1] = y;
		mMax[2] = z;
	}

	void include(const RmReal *v)
	{
		if ( v[0] < mMin[0] ) mMin[0] = v[0];
		if ( v[1] < mMin[1] ) mMin[1] = v[1];
		if ( v[2] < mMin[2] ) mMin[2] = v[2];

		if ( v[0] > mMax[0] ) mMax[0] = v[0];
		if ( v[1] > mMax[1] ) mMax[1] = v[1];
		if ( v[2] > mMax[2] ) mMax[2] = v[2];
	}

	void getCenter(RmReal *center) const
	{
		center[0] = (mMin[0]+mMax[0])*0.5f;
		center[1] = (mMin[1]+mMax[1])*0.5f;
		center[2] = (mMin[2]+mMax[2])*0.5f;
	}

	bool containsTriangle(const RmReal *p1,const RmReal *p2,const RmReal *p3) const
	{
		RmReal boxCenter[3];
		RmReal	boxHalfSize[3];
		RmReal triVerts[3][3];

		boxCenter[0] = (mMin[0]+mMax[0])*0.5f;
		boxCenter[1] = (mMin[1]+mMax[1])*0.5f;
		boxCenter[2] = (mMin[2]+mMax[2])*0.5f;

		boxHalfSize[0] = (mMax[0]-mMin[0])*0.5f;
		boxHalfSize[1] = (mMax[1]-mMin[1])*0.5f;
		boxHalfSize[2] = (mMax[2]-mMin[2])*0.5f;

		triVerts[0][0] = p1[0];
		triVerts[0][1] = p1[1];
		triVerts[0][2] = p1[2];

		triVerts[1][0] = p2[0];
		triVerts[1][1] = p2[1];
		triVerts[1][2] = p2[2];

		triVerts[2][0] = p3[0];
		triVerts[2][1] = p3[1];
		triVerts[2][2] = p3[2];

		int ret = triBoxOverlap(boxCenter,boxHalfSize,triVerts);

		return ret == 1 ? true : false;
	}
	RmUint32 clipTestXYZ(const RmReal *p) const
	{
		RmUint32 ocode = 0;
		if ( p[0] < mMin[0] ) ocode|=OLEFT;
		if ( p[0] > mMax[0] ) ocode|=ORIGHT;

		if ( p[1] < mMin[1] ) ocode|=OTOP;
		if ( p[1] > mMax[1] ) ocode|=OBOTTOM;

		if ( p[2] < mMin[2] ) ocode|=OFRONT;
		if ( p[2] > mMax[2] ) ocode|=OBACK;

		return ocode;
	};


	bool containsLineSegment(const RmReal *p1,const RmReal *p2,RmUint32 &acode) const
	{
		acode = 0;
		RmUint32 ocode1 = clipTestXYZ(p1);
		if ( !ocode1 ) return true;
		RmUint32 ocode2 = clipTestXYZ(p2);
		if ( !ocode2 ) return true;
		acode = ocode1 & ocode2;
		if ( acode ) return false;
		acode = acode;
		return true;
	}

	void clamp(const BoundsAABB &aabb)
	{
		if ( mMin[0] < aabb.mMin[0] ) mMin[0] = aabb.mMin[0];
		if ( mMin[1] < aabb.mMin[1] ) mMin[1] = aabb.mMin[1];
		if ( mMin[2] < aabb.mMin[2] ) mMin[2] = aabb.mMin[2];
		if ( mMax[0] > aabb.mMax[0] ) mMax[0] = aabb.mMax[0];
		if ( mMax[1] > aabb.mMax[1] ) mMax[1] = aabb.mMax[1];
		if ( mMax[2] > aabb.mMax[2] ) mMax[2] = aabb.mMax[2];
	}



	RmReal		mMin[3];
	RmReal		mMax[3];
};


class NodeAABB;

class NodeInterface
{
public:
	virtual NodeAABB * getNode(void) = 0;
	virtual void getFaceNormal(RmUint32 tri,RmReal *faceNormal) = 0;
};





	class NodeAABB
	{
	public:
		NodeAABB(void)
		{
			mLeft = NULL;
			mRight = NULL;
			mLeafTriangleIndex= TRI_EOF;
		}

		NodeAABB(RmUint32 vcount,const RmReal *vertices,RmUint32 tcount,RmUint32 *indices,
			RmUint32 maxDepth,	// Maximum recursion depth for the triangle mesh.
			RmUint32 minLeafSize,	// minimum triangles to treat as a 'leaf' node.
			RmReal	minAxisSize,
			NodeInterface *callback,
			TriVector &leafTriangles)	// once a particular axis is less than this size, stop sub-dividing.

		{
			mLeft = NULL;
			mRight = NULL;
			mLeafTriangleIndex = TRI_EOF;
			TriVector triangles;
			triangles.reserve(tcount);
			for (RmUint32 i=0; i<tcount; i++)
			{
				triangles.push_back(i);
			}
			mBounds.setMin( vertices );
			mBounds.setMax( vertices );
			const RmReal *vtx = vertices+3;
			for (RmUint32 i=1; i<vcount; i++)
			{
				mBounds.include( vtx );
				vtx+=3;
			}
			split(triangles,vcount,vertices,tcount,indices,0,maxDepth,minLeafSize,minAxisSize,callback,leafTriangles);
		}

		NodeAABB(const BoundsAABB &aabb)
		{
			mBounds = aabb;
			mLeft = NULL;
			mRight = NULL;
			mLeafTriangleIndex = TRI_EOF;
		}

		~NodeAABB(void)
		{
		}

		// here is where we split the mesh..
		void split(const TriVector &triangles,
			RmUint32 vcount,
			const RmReal *vertices,
			RmUint32 tcount,
			const RmUint32 *indices,
			RmUint32 depth,
			RmUint32 maxDepth,	// Maximum recursion depth for the triangle mesh.
			RmUint32 minLeafSize,	// minimum triangles to treat as a 'leaf' node.
			RmReal	minAxisSize,
			NodeInterface *callback,
			TriVector &leafTriangles)	// once a particular axis is less than this size, stop sub-dividing.

		{
			// Find the longest axis of the bounding volume of this node
			RmReal dx = mBounds.mMax[0] - mBounds.mMin[0];
			RmReal dy = mBounds.mMax[1] - mBounds.mMin[1];
			RmReal dz = mBounds.mMax[2] - mBounds.mMin[2];

			AxisAABB axis = AABB_XAXIS;
			RmReal laxis = dx;
			if ( dy > dx )
			{
				axis = AABB_YAXIS;
				laxis = dy;
			}
			if ( dz > dx && dz > dy )
			{
				axis = AABB_ZAXIS;
				laxis = dz;
			}

			RmUint32 count = triangles.size();

			// if the number of triangles is less than the minimum allowed for a leaf node or...
			// we have reached the maximum recursion depth or..
			// the width of the longest axis is less than the minimum axis size then...
			// we create the leaf node and copy the triangles into the leaf node triangle array.
			if ( count < minLeafSize || depth >= maxDepth || laxis < minAxisSize )
			{ 
				// Copy the triangle indices into the leaf triangles array
				mLeafTriangleIndex = leafTriangles.size(); // assign the array start location for these leaf triangles.
				leafTriangles.push_back(count);
				for (TriVector::const_iterator i=triangles.begin(); i!=triangles.end(); ++i)
				{
					leafTriangles.push_back( *i );
				}
			}
			else
			{
				RmReal center[3];
				mBounds.getCenter(center);
				BoundsAABB b1,b2;
				splitRect(axis,mBounds,b1,b2,center);

				// Compute two bounding boxes based upon the split of the longest axis

				BoundsAABB leftBounds,rightBounds;

				TriVector leftTriangles;
				TriVector rightTriangles;


				// Create two arrays; one of all triangles which intersect the 'left' half of the bounding volume node
				// and another array that includes all triangles which intersect the 'right' half of the bounding volume node.
				for (TriVector::const_iterator i=triangles.begin(); i!=triangles.end(); ++i)
				{
					RmUint32 tri = (*i); 

					{
						RmUint32 i1 = indices[tri*3+0];
						RmUint32 i2 = indices[tri*3+1];
						RmUint32 i3 = indices[tri*3+2];

						const RmReal *p1 = &vertices[i1*3];
						const RmReal *p2 = &vertices[i2*3];
						const RmReal *p3 = &vertices[i3*3];

						if ( b1.containsTriangle(p1,p2,p3))
						{
							if ( leftTriangles.empty() )
							{
								leftBounds.setMin(p1);
								leftBounds.setMax(p1);
							}
							leftBounds.include(p1);
							leftBounds.include(p2);
							leftBounds.include(p3);
							leftTriangles.push_back(tri); // Add this triangle to the 'left triangles' array and revise the left triangles bounding volume
						}

						if ( b2.containsTriangle(p1,p2,p3))
						{
							if ( rightTriangles.empty() )
							{
								rightBounds.setMin(p1);
								rightBounds.setMax(p1);
							}
							rightBounds.include(p1);
							rightBounds.include(p2);
							rightBounds.include(p3);
							rightTriangles.push_back(tri); // Add this triangle to the 'right triangles' array and revise the right triangles bounding volume.
						}
					}
				}

				if ( !leftTriangles.empty() ) // If there are triangles in the left half then...
				{
					leftBounds.clamp(b1); // we have to clamp the bounding volume so it stays inside the parent volume.
					mLeft = callback->getNode();	// get a new AABB node
					new ( mLeft ) NodeAABB(leftBounds);		// initialize it to default constructor values.  
					// Then recursively split this node.
					mLeft->split(leftTriangles,vcount,vertices,tcount,indices,depth+1,maxDepth,minLeafSize,minAxisSize,callback,leafTriangles);
				}

				if ( !rightTriangles.empty() ) // If there are triangles in the right half then..
				{
					rightBounds.clamp(b2);	// clamps the bounding volume so it stays restricted to the size of the parent volume.
					mRight = callback->getNode(); // allocate and default initialize a new node
					new ( mRight ) NodeAABB(rightBounds);
					// Recursively split this node.
					mRight->split(rightTriangles,vcount,vertices,tcount,indices,depth+1,maxDepth,minLeafSize,minAxisSize,callback,leafTriangles);
				}

			}
		}

		void splitRect(AxisAABB axis,const BoundsAABB &source,BoundsAABB &b1,BoundsAABB &b2,const RmReal *midpoint)
		{
			switch ( axis )
			{
				case AABB_XAXIS:
					{
						b1.setMin( source.mMin );
						b1.setMax( midpoint[0], source.mMax[1], source.mMax[2] );

						b2.setMin( midpoint[0], source.mMin[1], source.mMin[2] );
						b2.setMax(source.mMax);
					}
					break;
				case AABB_YAXIS:
					{
						b1.setMin(source.mMin);
						b1.setMax(source.mMax[0], midpoint[1], source.mMax[2]);

						b2.setMin(source.mMin[0], midpoint[1], source.mMin[2]);
						b2.setMax(source.mMax);
					}
					break;
				case AABB_ZAXIS:
					{
						b1.setMin(source.mMin);
						b1.setMax(source.mMax[0], source.mMax[1], midpoint[2]);

						b2.setMin(source.mMin[0], source.mMin[1], midpoint[2]);
						b2.setMax(source.mMax);
					}
					break;
			}
		}


		virtual void raycast(bool &hit,
							const RmReal *from,
							const RmReal *to,
							const RmReal *dir,
							RmReal *hitLocation,
							RmReal *hitNormal,
							RmReal *hitDistance,
							const RmReal *vertices,
							const RmUint32 *indices,
							RmReal &nearestDistance,
							NodeInterface *callback,
							unsigned char *raycastTriangles,
							unsigned char raycastFrame,
							const TriVector &leafTriangles)
		{
			unsigned acode;
			if ( !mBounds.containsLineSegment(from,to,acode))
			{
				return;
			}
			if ( acode )
			{
				RmReal sect[3];
				RmReal nd = nearestDistance;
				if ( !intersectLineSegmentAABB(mBounds.mMin,mBounds.mMax,from,dir,nd,sect) )
				{
					return;	
				}
			}

			if ( mLeafTriangleIndex != TRI_EOF )
			{
				const RmUint32 *scan = &leafTriangles[mLeafTriangleIndex];
				RmUint32 count = *scan++;
				for (RmUint32 i=0; i<count; i++)
				{
					RmUint32 tri = *scan++;
					if ( raycastTriangles[tri] != raycastFrame )
					{
						raycastTriangles[tri] = raycastFrame;
						RmUint32 i1 = indices[tri*3+0];
						RmUint32 i2 = indices[tri*3+1];
						RmUint32 i3 = indices[tri*3+2];

						const RmReal *p1 = &vertices[i1*3];
						const RmReal *p2 = &vertices[i2*3];
						const RmReal *p3 = &vertices[i3*3];

						RmReal t;
						if ( rayIntersectsTriangle(from,dir,p1,p2,p3,t))
						{
							if ( t < nearestDistance )
							{
								nearestDistance = t;
								if ( hitLocation )
								{
									hitLocation[0] = from[0]+dir[0]*t;
									hitLocation[1] = from[1]+dir[1]*t;
									hitLocation[2] = from[2]+dir[2]*t;
								}
								if ( hitNormal )
								{
									callback->getFaceNormal(tri,hitNormal);
								}
								if ( hitDistance )
								{
									*hitDistance = t;
								}
								hit = true;
							}
						}
					}
				}
			}
			else
			{
				if ( mLeft )
				{
					mLeft->raycast(hit,from,to,dir,hitLocation,hitNormal,hitDistance,vertices,indices,nearestDistance,callback,raycastTriangles,raycastFrame,leafTriangles);
				}
				if ( mRight )
				{
					mRight->raycast(hit,from,to,dir,hitLocation,hitNormal,hitDistance,vertices,indices,nearestDistance,callback,raycastTriangles,raycastFrame,leafTriangles);
				}
			}
		}

		NodeAABB		*mLeft;			// left node
		NodeAABB		*mRight;		// right node
		BoundsAABB		mBounds;		// bounding volume of node
		RmUint32		mLeafTriangleIndex;	// if it is a leaf node; then these are the triangle indices.
	};

class MyRaycastMesh : public RaycastMesh, public NodeInterface
{
public:

	MyRaycastMesh(RmUint32 vcount,const RmReal *vertices,RmUint32 tcount,const RmUint32 *indices,RmUint32 maxDepth,RmUint32 minLeafSize,RmReal minAxisSize)
	{
		mRaycastFrame = 0;
		if ( maxDepth < 2 )
		{
			maxDepth = 2;
		}
		if ( maxDepth > 15 )
		{
			maxDepth = 15;
		}
		RmUint32 pow2Table[16] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 65536 };
		mMaxNodeCount = 0;
		for (RmUint32 i=0; i<=maxDepth; i++)
		{
			mMaxNodeCount+=pow2Table[i];
		}
		mNodes = new NodeAABB[mMaxNodeCount];
		mNodeCount = 0;
		mVcount = vcount;
		mVertices = (RmReal *)::malloc(sizeof(RmReal)*3*vcount);
		memcpy(mVertices,vertices,sizeof(RmReal)*3*vcount);
		mTcount = tcount;
		mIndices = (RmUint32 *)::malloc(sizeof(RmUint32)*tcount*3);
		memcpy(mIndices,indices,sizeof(RmUint32)*tcount*3);
		mRaycastTriangles = (unsigned char *)::malloc(tcount);
		memset(mRaycastTriangles,0,tcount);
		mRoot = getNode();
		mFaceNormals = NULL;
		new ( mRoot ) NodeAABB(mVcount,mVertices,mTcount,mIndices,maxDepth,minLeafSize,minAxisSize,this,mLeafTriangles);
	}

	~MyRaycastMesh(void)
	{
		delete []mNodes;
		::free(mVertices);
		::free(mIndices);
		::free(mFaceNormals);
		::free(mRaycastTriangles);
	}

	virtual bool raycast(const RmReal *from,const RmReal *to,RmReal *hitLocation,RmReal *hitNormal,RmReal *hitDistance)
	{
		bool ret = false;

		RmReal dir[3];
		dir[0] = to[0] - from[0];
		dir[1] = to[1] - from[1];
		dir[2] = to[2] - from[2];
		RmReal distance = sqrtf( dir[0]*dir[0] + dir[1]*dir[1]+dir[2]*dir[2] );
		if ( distance < 0.0000000001f ) return false;
		RmReal recipDistance = 1.0f / distance;
		dir[0]*=recipDistance;
		dir[1]*=recipDistance;
		dir[2]*=recipDistance;
#if USE_BRUTE_FORCE
		const RmUint32 *indices = mIndices;
		const RmReal *vertices = mVertices;
		RmReal nearestDistance = distance;
		for (RmUint32 tri=0; tri<mTcount; tri++)
		{
			RmUint32 i1 = indices[tri*3+0];
			RmUint32 i2 = indices[tri*3+1];
			RmUint32 i3 = indices[tri*3+2];

			const RmReal *p1 = &vertices[i1*3];
			const RmReal *p2 = &vertices[i2*3];
			const RmReal *p3 = &vertices[i3*3];

			RmReal t;
			if ( rayIntersectsTriangle(from,dir,p1,p2,p3,t))
			{
				if ( t < nearestDistance )
				{
					nearestDistance = t;
					if ( hitLocation )
					{
						hitLocation[0] = from[0]+dir[0]*t;
						hitLocation[1] = from[1]+dir[1]*t;
						hitLocation[2] = from[2]+dir[2]*t;
					}
					if ( hitNormal )
					{
						getFaceNormal(tri,hitNormal);
					}
					if ( hitDistance )
					{
						*hitDistance = t;
					}
					ret = true;
				}
			}
		}
#else
		mRaycastFrame++;
		mRoot->raycast(ret,from,to,dir,hitLocation,hitNormal,hitDistance,mVertices,mIndices,distance,this,mRaycastTriangles,mRaycastFrame,mLeafTriangles);
#endif


		return ret;
	}

	virtual void release(void)
	{
		delete this;
	}

	virtual const RmReal * getBoundMin(void) const // return the minimum bounding box
	{
		return mRoot->mBounds.mMin;
	}
	virtual const RmReal * getBoundMax(void) const // return the maximum bounding box.
	{
		return mRoot->mBounds.mMax;
	}

	virtual NodeAABB * getNode(void) 
	{
		assert( mNodeCount < mMaxNodeCount );
		NodeAABB *ret = &mNodes[mNodeCount];
		mNodeCount++;
		return ret;
	}

	virtual void getFaceNormal(RmUint32 tri,RmReal *faceNormal) 
	{
		if ( mFaceNormals == NULL )
		{
			mFaceNormals = (RmReal *)::malloc(sizeof(RmReal)*3*mTcount);
			for (RmUint32 i=0; i<mTcount; i++)
			{
				RmUint32 i1 = mIndices[i*3+0];
				RmUint32 i2 = mIndices[i*3+1];
				RmUint32 i3 = mIndices[i*3+2];
				const RmReal*p1 = &mVertices[i1*3];
				const RmReal*p2 = &mVertices[i2*3];
				const RmReal*p3 = &mVertices[i3*3];
				RmReal *dest = &mFaceNormals[i*3];
				computePlane(p3,p2,p1,dest);
			}
		}
		const RmReal *src = &mFaceNormals[tri*3];
		faceNormal[0] = src[0];
		faceNormal[1] = src[1];
		faceNormal[2] = src[2];
	}

	unsigned char	mRaycastFrame;
	unsigned char	*mRaycastTriangles;
	RmUint32		mVcount;
	RmReal			*mVertices;
	RmReal			*mFaceNormals;
	RmUint32		mTcount;
	RmUint32		*mIndices;
	NodeAABB		*mRoot;
	RmUint32		mNodeCount;
	RmUint32		mMaxNodeCount;
	NodeAABB		*mNodes;
	TriVector		mLeafTriangles;
};

};



using namespace RAYCAST_MESH;


RaycastMesh * createRaycastMesh(RmUint32 vcount,		// The number of vertices in the source triangle mesh
								const RmReal *vertices,		// The array of vertex positions in the format x1,y1,z1..x2,y2,z2.. etc.
								RmUint32 tcount,		// The number of triangles in the source triangle mesh
								const RmUint32 *indices, // The triangle indices in the format of i1,i2,i3 ... i4,i5,i6, ...
								RmUint32 maxDepth,	// Maximum recursion depth for the triangle mesh.
								RmUint32 minLeafSize,	// minimum triangles to treat as a 'leaf' node.
								RmReal	minAxisSize	// once a particular axis is less than this size, stop sub-dividing.
								)
{
	MyRaycastMesh *m = new MyRaycastMesh(vcount,vertices,tcount,indices,maxDepth,minLeafSize,minAxisSize);
	return static_cast< RaycastMesh * >(m);
}



