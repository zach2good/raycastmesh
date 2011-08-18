#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#pragma warning(disable:4996)

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "RaycastMesh.h"
#include "wavefront.h"

int main(int argc,const char **argv)
{
	if ( argc != 2 )
	{
		printf("Usage: RaycastMesh <wavefront.obj>\r\n");
	}
	else
	{
		printf("Loading wavefront obj file '%s'\r\n", argv[1] );
		WavefrontObj obj;
		unsigned int tcount = obj.loadObj(argv[1],false);
		if ( tcount )
		{
			RaycastMesh *rm = createRaycastMesh(obj.mVertexCount,obj.mVertices,obj.mTriCount,(const RmUint32 *)obj.mIndices);
			if ( rm )
			{
#define IMAGE_SIZE 1024
				const float *bmin = rm->getBoundMin();
				const float *bmax = rm->getBoundMax();
				unsigned char *image = (unsigned char *)::malloc(IMAGE_SIZE*IMAGE_SIZE);
				RmReal *distances = (RmReal *)::malloc(IMAGE_SIZE*IMAGE_SIZE*sizeof(RmReal));
				printf("Performing %d raycasts.\r\n", IMAGE_SIZE*IMAGE_SIZE);
				RmReal minDistance = 1e9;
				RmReal maxDistance = -1e9;
				for (RmUint32 y=0; y<IMAGE_SIZE; y++)
				{
					float fractionY = (float)(IMAGE_SIZE-y)/IMAGE_SIZE;
					float fy = (bmax[1]-bmin[1])*fractionY+bmin[1];
					for (RmUint32 x=0; x<IMAGE_SIZE; x++)
					{
						float fractionZ = (float)x/IMAGE_SIZE;
						float fz = (bmax[2]-bmin[2])*fractionZ+bmin[2];

						RmReal from[3] = { -1000, fy, fz };
						RmReal to[3] = { 1000, fy, fz };
						RmReal hitLocation[3];
						RmReal normal[3];
						RmReal hitDistance;
						bool hit = rm->raycast(from,to,hitLocation,normal,&hitDistance);
						if ( hit )
						{
							if ( hitDistance < minDistance ) minDistance = hitDistance;
							if ( hitDistance > maxDistance ) maxDistance = hitDistance;
							distances[y*IMAGE_SIZE+x] = hitDistance;
							image[y*IMAGE_SIZE+x] = 0;
						}
						else
						{
							image[y*IMAGE_SIZE+x] = 255;
						}
					}
				}

				printf("Finished raycasting; converting image to grayscale.\r\n");
				RmReal grayScale = 255.0f / (maxDistance-minDistance);
				for (RmUint32 y=0; y<IMAGE_SIZE; y++)
				{
					for (RmUint32 x=0; x<IMAGE_SIZE; x++)
					{
						if ( image[y*IMAGE_SIZE+x] == 0 )
						{
							RmReal d = distances[y*IMAGE_SIZE+x]-minDistance;
							d = 255-(d*grayScale);
							image[y*IMAGE_SIZE+x] = (unsigned char)d;
						}
					}
				}

				printf("Saving image file RaycastMesh.png\r\n");
				stbi_write_png("RaycastMesh.png",IMAGE_SIZE,IMAGE_SIZE,1,image,IMAGE_SIZE);

				::free(image);
				::free(distances);

				rm->release();
			}

		}
		else
		{
			printf("Failed to find any triangles in the input mesh.\r\n");
		}
	}
}
