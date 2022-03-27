#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <numeric>
#include <vector>
#include <algorithm>

#pragma warning(disable : 4996)

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "RaycastMesh.h"
#include "stb_image_write.h"
#include "wavefront.h"

#define BRUTE_FORCE_VALIDATION 1 // enable this to do exact validation of the AABB raycasting system
#define RAYCAST_TOP_DOWN       0

using namespace std::chrono;

int main(int argc, const char** argv)
{
    if (argc != 2)
    {
        printf("Usage: RaycastMesh <wavefront.obj>\n");
    }
    else
    {
        std::vector<unsigned long long> seconds(2000000);

        printf("Loading wavefront obj file '%s'\n", argv[1]);
        WavefrontObj obj;
        unsigned int tcount = obj.loadObj(argv[1], false);
        if (tcount)
        {
            printf("Creating Raycast Mesh. Building AABB\n");
            RaycastMesh* rm = createRaycastMesh(obj.mVertexCount, obj.mVertices, obj.mTriCount, (const RmUint32*)obj.mIndices);
            if (rm)
            {
#define IMAGE_SIZE 256

                const float* bmin = rm->getBoundMin();
                const float* bmax = rm->getBoundMax();

                unsigned char* image     = (unsigned char*)::malloc(IMAGE_SIZE * IMAGE_SIZE);
                RmReal*        distances = (RmReal*)::malloc(IMAGE_SIZE * IMAGE_SIZE * sizeof(RmReal));
                printf("Performing %d raycasts.\n", IMAGE_SIZE * IMAGE_SIZE);
                RmReal minDistance = 1e9;
                RmReal maxDistance = -1e9;
                for (RmUint32 y = 0; y < IMAGE_SIZE; y++)
                {
#if RAYCAST_TOP_DOWN
                    float fractionZ = (float)(IMAGE_SIZE - y) / IMAGE_SIZE;
                    float fz        = (bmax[2] - bmin[2]) * fractionZ + bmin[2];
#else
                    float fractionY = (float)(IMAGE_SIZE - y) / IMAGE_SIZE;
                    float fy        = (bmax[1] - bmin[1]) * fractionY + bmin[1];
#endif
                    for (RmUint32 x = 0; x < IMAGE_SIZE; x++)
                    {
#if RAYCAST_TOP_DOWN
                        float  fractionX = (float)x / IMAGE_SIZE;
                        float  fx        = (bmax[0] - bmin[0]) * fractionX + bmin[0];
                        RmReal from[3]   = { fx, 1000, fz };
                        RmReal to[3]     = { fx, -1000, fz };
#else
                        float  fractionZ = (float)x / IMAGE_SIZE;
                        float  fz        = (bmax[2] - bmin[2]) * fractionZ + bmin[2];
                        RmReal from[3]   = { -1000, fy, fz };
                        RmReal to[3]     = { 1000, fy, fz };
#endif
                        RmReal hitLocation[3];
                        RmReal normal[3];
                        RmReal hitDistance;

                        auto start = high_resolution_clock::now();

                        bool hit = rm->raycast(from, to, hitLocation, normal, &hitDistance);

                        auto stop = high_resolution_clock::now();

                        seconds.emplace_back(duration_cast<nanoseconds>(stop - start).count());

// If 'BRUTE_FORCE_VALIDATION' is enabled then we manually check every single triangle and make sure we get the same result as the optimized code path
#if BRUTE_FORCE_VALIDATION
                        {
                            RmReal bruteHitLocation[3];
                            RmReal bruteNormal[3];
                            RmReal bruteHitDistance;
                            bool   bruteHit = rm->bruteForceRaycast(from, to, bruteHitLocation, bruteNormal, &bruteHitDistance);
                            assert(bruteHit == hit);
                            if (hit)
                            {
                                assert(bruteHitDistance == hitDistance);
                                assert(bruteNormal[0] == normal[0]);
                                assert(bruteNormal[1] == normal[1]);
                                assert(bruteNormal[2] == normal[2]);
                                assert(bruteHitLocation[0] == hitLocation[0]);
                                assert(bruteHitLocation[1] == hitLocation[1]);
                                assert(bruteHitLocation[2] == hitLocation[2]);
                            }
                        }
#endif
                        if (hit)
                        {
                            if (hitDistance < minDistance)
                            {
                                minDistance = hitDistance;
                            }
                            if (hitDistance > maxDistance)
                            {
                                maxDistance = hitDistance;
                            }
                            distances[y * IMAGE_SIZE + x] = hitDistance;
                            image[y * IMAGE_SIZE + x]     = 0;
                        }
                        else
                        {
                            image[y * IMAGE_SIZE + x] = 255;
                        }
                    }
                }

                printf("Finished raycasting; converting image to grayscale.\n");
                RmReal grayScale = 255.0f / (maxDistance - minDistance);
                for (RmUint32 y = 0; y < IMAGE_SIZE; y++)
                {
                    for (RmUint32 x = 0; x < IMAGE_SIZE; x++)
                    {
                        if (image[y * IMAGE_SIZE + x] == 0)
                        {
                            RmReal d                  = distances[y * IMAGE_SIZE + x] - minDistance;
                            d                         = 255 - (d * grayScale);
                            image[y * IMAGE_SIZE + x] = (unsigned char)d;
                        }
                    }
                }

                printf("Saving image file RaycastMesh.png\n");
                stbi_write_png("RaycastMesh.png", IMAGE_SIZE, IMAGE_SIZE, 1, image, IMAGE_SIZE);

                printf("Average Time: %lli ns\n", std::reduce(seconds.begin(), seconds.end()) / seconds.size());
				auto mixmax = std::minmax_element(seconds.begin(), seconds.end());
				printf("Min Time: %lli ns\n", *mixmax.first);
				printf("Max Time: %lli ns\n", *mixmax.second);

                ::free(image);
                ::free(distances);

                rm->release();
            }
        }
        else
        {
            printf("Failed to find any triangles in the input mesh.\n");
        }
    }
}
