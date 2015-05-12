
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// accelerators/bvh.cpp*
#include "stdafx.h"
#include "accelerators/bvh.h"
#include "probes.h"
#include "paramset.h"
#define FLT_MAX 3.40282347E+38F
#define MORTONBIT 10
#define MAXDIGIT 10
#define BASE 8

// BVHAccel Local Declarations
struct BVHPrimitiveInfo {
    BVHPrimitiveInfo() { }
    BVHPrimitiveInfo(int pn, const BBox &b)
        : primitiveNumber(pn), bounds(b) {
        centroid = .5f * b.pMin + .5f * b.pMax;
    }
    int primitiveNumber;
    Point centroid;
    BBox bounds;
    uint32_t mortonCode;
};

struct MortonSort {
    inline bool operator () (const BVHPrimitiveInfo& p1, const BVHPrimitiveInfo& p2) {
        return ( p1.mortonCode < p2.mortonCode );
    }
};

struct BVHBuildNode {
    // BVHBuildNode Public Methods
    BVHBuildNode() { children[0] = children[1] = NULL; }
    void InitLeaf(uint32_t first, uint32_t n, const BBox &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
    }
    void InitInterior(uint32_t axis, BVHBuildNode *c0, BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
    }
    BBox bounds;
    BVHBuildNode *children[2];
    uint32_t splitAxis, firstPrimOffset, nPrimitives;
};


struct CompareToMid {
    CompareToMid(int d, float m) { dim = d; mid = m; }
    int dim;
    float mid;
    bool operator()(const BVHPrimitiveInfo &a) const {
        return a.centroid[dim] < mid;
    }
};


struct ComparePoints {
    ComparePoints(int d) { dim = d; }
    int dim;
    bool operator()(const BVHPrimitiveInfo &a,
                    const BVHPrimitiveInfo &b) const {
        return a.centroid[dim] < b.centroid[dim];
    }
};


struct CompareToBucket {
    CompareToBucket(int split, int num, int d, const BBox &b)
        : centroidBounds(b)
    { splitBucket = split; nBuckets = num; dim = d; }
    bool operator()(const BVHPrimitiveInfo &p) const;

    int splitBucket, nBuckets, dim;
    const BBox &centroidBounds;
};


bool CompareToBucket::operator()(const BVHPrimitiveInfo &p) const {
    int b = nBuckets * ((p.centroid[dim] - centroidBounds.pMin[dim]) /
            (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
    if (b == nBuckets) b = nBuckets-1;
    Assert(b >= 0 && b < nBuckets);
    return b <= splitBucket;
}


struct LinearBVHNode {
    BBox bounds;
    union {
        uint32_t primitivesOffset;    // leaf
        uint32_t secondChildOffset;   // interior
    };

    uint8_t nPrimitives;  // 0 -> interior node
    uint8_t axis;         // interior node: xyz
    uint8_t pad[2];       // ensure 32 byte total size
};

static inline bool IntersectP(const BBox &bounds, const Ray &ray,
        const Vector &invDir, const uint32_t dirIsNeg[3]) {
    // Check for ray intersection against $x$ and $y$ slabs
    float tmin =  (bounds[  dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tmax =  (bounds[1-dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tymin = (bounds[  dirIsNeg[1]].y - ray.o.y) * invDir.y;
    float tymax = (bounds[1-dirIsNeg[1]].y - ray.o.y) * invDir.y;
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    // Check for ray intersection against $z$ slab
    float tzmin = (bounds[  dirIsNeg[2]].z - ray.o.z) * invDir.z;
    float tzmax = (bounds[1-dirIsNeg[2]].z - ray.o.z) * invDir.z;
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    return (tmin < ray.maxt) && (tmax > ray.mint);
}

uint32_t part1By2(uint32_t x) {
    // REFERENCE: https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
    x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
    x = (x ^ (x << 16)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x ^ (x <<  8)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x ^ (x <<  4)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x ^ (x <<  2)) & 0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    return x;
}

uint32_t encodeMorton3(uint32_t x, uint32_t y, uint32_t z) {
    // REFERENCE: https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
    return ( part1By2(z) << 2 ) + ( part1By2(y) << 1 ) + part1By2( x );
}
static void computeMortonCodes(vector<BVHPrimitiveInfo>& buildData) {
    if ( buildData.size() == 0 ) return;
    // compute the bounding box for all primitive centroids;
    size_t numPrimitives = buildData.size();
    BBox centroidBounds;
    for (uint32_t i = 0; i < numPrimitives; ++i)
        centroidBounds = Union(centroidBounds, buildData[i].centroid);
    
    // compute relative coordinates and mortoncodes
    for (uint32_t i = 0; i < numPrimitives; i ++) {
        uint32_t rel_x = (uint32_t)(( buildData[i].centroid.x - centroidBounds.pMin.x ) /
                                    ( centroidBounds.pMax.x - centroidBounds.pMin.x ) * (1 << MORTONBIT));
        uint32_t rel_y = (uint32_t)(( buildData[i].centroid.y - centroidBounds.pMin.y ) /
                                    ( centroidBounds.pMax.y - centroidBounds.pMin.y ) * (1 << MORTONBIT));
        uint32_t rel_z = (uint32_t)(( buildData[i].centroid.z - centroidBounds.pMin.z ) /
                                    ( centroidBounds.pMax.z - centroidBounds.pMin.z ) * (1 << MORTONBIT));
        rel_x -= ( rel_x == (1 << MORTONBIT));
        rel_y -= ( rel_y == (1 << MORTONBIT));
        rel_z -= ( rel_z == (1 << MORTONBIT));
        buildData[i].mortonCode = encodeMorton3(rel_x, rel_y, rel_z);
    }
}

static void countSort( vector<BVHPrimitiveInfo>& arr, int exp ) {
    // sort number from 0 to range
    int count[BASE] = {0};
    vector<BVHPrimitiveInfo> output(arr.size());
    for (int i = 0; i < arr.size(); i++) {
        count[(arr[i].mortonCode / exp) % BASE] ++;
    }
    for (int i = 1; i < BASE; i++) {
        count[i] += count[i - 1];
    }
    for (int i = arr.size() - 1; i >= 0; i--) {
        output[count[(arr[i].mortonCode / exp) % BASE] - 1] = arr[i];
        count[(arr[i].mortonCode / exp) % BASE] --;
    }
    arr = output;
}

static void radixSort( vector<BVHPrimitiveInfo>& arr ) {
    int exponent = 1;
    for (int d = 0; d < MAXDIGIT; d++) {
        countSort(arr, exponent);
        exponent *= BASE;
    }
}

static void sortPrimitivesByMortonCodes(vector<BVHPrimitiveInfo>& buildData) {
    //radixSort(buildData);
    std::sort(buildData.begin(), buildData.end(), MortonSort());
}

inline uint32_t countReductFunc(uint32_t nPrimitives) {
    if ((float_t) nPrimitives < 0) {
        printf("Num of primitives are too large \n");
        return 0;
    }
    return 0.5 * pow(4, 0.5 + 0.2) * pow(nPrimitives, 0.5 - 0.2) + 0.5;
}

inline float_t getDist(BVHBuildNode *node1, BVHBuildNode* node2) {
    BBox box = Union(node1->bounds, node2->bounds);
    return box.SurfaceArea();
}

struct ClusterDist {
    uint32_t bestMatchIdx;
    float_t minDist;
    ClusterDist(uint32_t idx, float_t d): bestMatchIdx(idx), minDist(d) {}
};

static ClusterDist findBestMatch(const vector<BVHBuildNode *>& clusterList, uint32_t clusterIdx ) {
    //Assert(clusterIdx < clusterList.size()); // TODO: delete
    if (clusterList.size() == 1) return  ClusterDist(clusterIdx, 0);
    float_t currDist, minDist = FLT_MAX;
    int bestMatch = -1;
    for (uint32_t i = 0; i < clusterList.size(); i++) {
        if (i == clusterIdx) continue;
        currDist = getDist(clusterList[i], clusterList[clusterIdx]);
        if (currDist < minDist) {
            minDist = currDist;
            bestMatch = i;
        }
    }
    return ClusterDist(bestMatch, minDist);
}

void BVHAccel::combineClusters(MemoryArena &buildArena, vector<BVHBuildNode *>& clusterList, uint32_t maxNumClusters, uint32_t *totalNodes) {
    if (clusterList.size() <= maxNumClusters ) return;
    
    // compute closest cluster index
    vector<ClusterDist> bestMatches;
    for (uint32_t i = 0; i < clusterList.size(); i++) {
        bestMatches.push_back(findBestMatch(clusterList, i));
    }
    
    // merge nodes
    int leftIdx = -1, rightIdx = -1;
    while (clusterList.size() > maxNumClusters) {
        // TODO: delete
//        printf("size of clusterList is %lu\n",clusterList.size());
//        for (int t = 0; t < clusterList.size(); t++) {
//            printf("min idx is %d \n", bestMatches[t].bestMatchIdx);
//        }
        float_t currDist, minDist = FLT_MAX;
        for (uint32_t i = 0; i < clusterList.size(); i++) {
            currDist = bestMatches[i].minDist;
            if (currDist < minDist) {
                minDist = currDist;
                //Assert(i < bestMatches[i].bestMatchIdx); // TODO: delete
                leftIdx = i;
                rightIdx = bestMatches[i].bestMatchIdx;
            }
        }
        //Assert(leftIdx >= 0 && rightIdx >= 0 && leftIdx < rightIdx);
        // merge two clusters
        BVHBuildNode *mergedCluster = buildArena.Alloc<BVHBuildNode>();
        mergedCluster->InitInterior(0, clusterList[leftIdx], clusterList[rightIdx]);
        (*totalNodes)++;
        
        // update cluster list
        clusterList[leftIdx] = mergedCluster;
        clusterList[rightIdx] = clusterList.back();
        clusterList.pop_back();

        // update closest cluster index
        bestMatches[leftIdx] = findBestMatch(clusterList, leftIdx);
        bestMatches[rightIdx] = bestMatches.back();
        bestMatches.pop_back();
        
        for (uint32_t j = 0; j < clusterList.size(); j++) {
            if (bestMatches[j].bestMatchIdx == leftIdx || bestMatches[j].bestMatchIdx == rightIdx || bestMatches[j].bestMatchIdx == bestMatches.size()) {
                bestMatches[j] = findBestMatch(clusterList, j);
            }
        }
    }
}

uint32_t makePartition(vector<BVHPrimitiveInfo> &buildData, uint32_t start, uint32_t end,  int32_t& curr_bit) {
    uint32_t s = start, t = end, mid;
    uint32_t mask = 1 << curr_bit;

    while (curr_bit >= 0) {
        if (!(buildData[s].mortonCode & mask) && (buildData[t].mortonCode & mask))
            break;
        curr_bit --;
        mask <<= 1;
    }
    if (curr_bit < 0)
        return (s + t) / 2;
    else {
        while (s < t) {
            mid = (s + t) / 2;
            if (buildData[mid].mortonCode & mask)
                t = mid;
            else
                s = mid + 1;
        }
    }
    // TODO: test
    return s;
}

vector<BVHBuildNode *> BVHAccel:: recursiveBuildAAC(MemoryArena &buildArena, vector<BVHPrimitiveInfo> &buildData, uint32_t start, uint32_t end, uint32_t *totalNodes, vector<Reference<Primitive> > &orderedPrims, int32_t curr_bit) {
    
    Assert(start != end);
    vector<BVHBuildNode *> clusterList;
    
    // base case
    uint32_t nPrimitives = end - start;
    if (nPrimitives < DeltaACC) {
        for (uint32_t i = start; i < end; ++i) {
            uint32_t firstPrimOffset = orderedPrims.size(); // TODO
            BVHBuildNode *leaf = buildArena.Alloc<BVHBuildNode>();
            (*totalNodes)++;
            leaf->InitLeaf(firstPrimOffset, 1, buildData[i].bounds);
            clusterList.push_back(leaf);
            uint32_t primNum = buildData[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        combineClusters(buildArena, clusterList, countReductFunc(DeltaACC), totalNodes);
        return clusterList;
    }
    else {
        uint32_t partitionPos = makePartition(buildData, start, end, curr_bit);
        vector<BVHBuildNode *> leftClusters = recursiveBuildAAC(buildArena, buildData, start, partitionPos, totalNodes, orderedPrims, curr_bit - 1);
        vector<BVHBuildNode *> rightClusters = recursiveBuildAAC(buildArena, buildData, partitionPos, end, totalNodes, orderedPrims, curr_bit - 1);
        
        leftClusters.reserve(leftClusters.size() + rightClusters.size());
        leftClusters.insert(leftClusters.end(), rightClusters.begin(), rightClusters.end());
        rightClusters.clear();
        combineClusters(buildArena, leftClusters, countReductFunc(nPrimitives), totalNodes);
        return leftClusters;
    }
}

BVHBuildNode *BVHAccel::buildAAC(MemoryArena &buildArena, vector<BVHPrimitiveInfo> &buildData, uint32_t start, uint32_t end, uint32_t *totalNodes, vector<Reference<Primitive> > &orderedPrims) {
    //Assert(start != end);
    //printf("There are %lu primitives.\n", buildData.size());
    // compute Morton Code for centers of P
    computeMortonCodes(buildData);
    sortPrimitivesByMortonCodes(buildData);

    // build the tree
    uint32_t curr_bit = 30;
    vector<BVHBuildNode *> clusters = recursiveBuildAAC(buildArena, buildData, start, end, totalNodes, orderedPrims, curr_bit);
    combineClusters(buildArena, clusters, 1, totalNodes);
    return clusters[0];
}


// BVHAccel Method Definitions
BVHAccel::BVHAccel(const vector<Reference<Primitive> > &p,
                   uint32_t mp, const string &sm) {
    maxPrimsInNode = min(255u, mp);
    for (uint32_t i = 0; i < p.size(); ++i)
        p[i]->FullyRefine(primitives);
    if (sm == "sah")         splitMethod = SPLIT_SAH;
    else if (sm == "middle") splitMethod = SPLIT_MIDDLE;
    else if (sm == "equal")  splitMethod = SPLIT_EQUAL_COUNTS;
    else if (sm == "aac") splitMethod = SPLIT_AAC;
    else {
        Warning("BVH split method \"%s\" unknown.  Using \"sah\".",
                sm.c_str());
        splitMethod = SPLIT_SAH;
    }

    if (primitives.size() == 0) {
        nodes = NULL;
        return;
    }
    // Build BVH from _primitives_
    PBRT_BVH_STARTED_CONSTRUCTION(this, primitives.size());

    // Initialize _buildData_ array for primitives
    vector<BVHPrimitiveInfo> buildData;
    buildData.reserve(primitives.size());
    for (uint32_t i = 0; i < primitives.size(); ++i) {
        BBox bbox = primitives[i]->WorldBound();
        buildData.push_back(BVHPrimitiveInfo(i, bbox));
    }

    // Recursively build BVH tree for primitives
    MemoryArena buildArena;
    uint32_t totalNodes = 0;
    vector<Reference<Primitive> > orderedPrims;
    orderedPrims.reserve(primitives.size());
    BVHBuildNode *root;
    if ( splitMethod == SPLIT_AAC) {
        root = buildAAC(buildArena, buildData, 0, primitives.size(), &totalNodes, orderedPrims);
    }
    else {
        root = recursiveBuild(buildArena, buildData, 0, primitives.size(),
                              &totalNodes, orderedPrims);
    }
    primitives.swap(orderedPrims);
        Info("BVH created with %d nodes for %d primitives (%.2f MB)", totalNodes,
             (int)primitives.size(), float(totalNodes * sizeof(LinearBVHNode))/(1024.f*1024.f));

    // Compute representation of depth-first traversal of BVH tree
    nodes = AllocAligned<LinearBVHNode>(totalNodes);
    for (uint32_t i = 0; i < totalNodes; ++i)
        new (&nodes[i]) LinearBVHNode;
    uint32_t offset = 0;
    flattenBVHTree(root, &offset);
    Assert(offset == totalNodes);
    PBRT_BVH_FINISHED_CONSTRUCTION(this);
}


BBox BVHAccel::WorldBound() const {
    return nodes ? nodes[0].bounds : BBox();
}


BVHBuildNode *BVHAccel::recursiveBuild(MemoryArena &buildArena,
        vector<BVHPrimitiveInfo> &buildData, uint32_t start,
        uint32_t end, uint32_t *totalNodes,
        vector<Reference<Primitive> > &orderedPrims) {
    Assert(start != end);
    (*totalNodes)++;
    BVHBuildNode *node = buildArena.Alloc<BVHBuildNode>();
    // Compute bounds of all primitives in BVH node
    BBox bbox;
    for (uint32_t i = start; i < end; ++i)
        bbox = Union(bbox, buildData[i].bounds);
    uint32_t nPrimitives = end - start;
    if (nPrimitives == 1) {
        // Create leaf _BVHBuildNode_
        uint32_t firstPrimOffset = orderedPrims.size();
        for (uint32_t i = start; i < end; ++i) {
            uint32_t primNum = buildData[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
    }
    else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        BBox centroidBounds;
        for (uint32_t i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, buildData[i].centroid);
        int dim = centroidBounds.MaximumExtent();

        // Partition primitives into two sets and build children
        uint32_t mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // If nPrimitives is no greater than maxPrimsInNode,
            // then all the nodes can be stored in a compact bvh node.
            if (nPrimitives <= maxPrimsInNode) {
                // Create leaf _BVHBuildNode_
                uint32_t firstPrimOffset = orderedPrims.size();
                for (uint32_t i = start; i < end; ++i) {
                    uint32_t primNum = buildData[i].primitiveNumber;
                    orderedPrims.push_back(primitives[primNum]);
                }
                node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
                return node;
            }
            else {
                // else if nPrimitives is greater than maxPrimsInNode, we
                // need to split it further to guarantee each node contains
                // no more than maxPrimsInNode primitives.
                node->InitInterior(dim,
                                   recursiveBuild(buildArena, buildData, start, mid,
                                                  totalNodes, orderedPrims),
                                   recursiveBuild(buildArena, buildData, mid, end,
                                                  totalNodes, orderedPrims));
                return node;
            }
        }

        // Partition primitives based on _splitMethod_
        switch (splitMethod) {
        case SPLIT_AAC: { break;} // TAG: ADDED
        case SPLIT_MIDDLE: {
            // Partition primitives through node's midpoint
            float pmid = .5f * (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]);
            BVHPrimitiveInfo *midPtr = std::partition(&buildData[start],
                                                      &buildData[end-1]+1,
                                                      CompareToMid(dim, pmid));
            mid = midPtr - &buildData[0];
            if (mid != start && mid != end)
                // for lots of prims with large overlapping bounding boxes, this
                // may fail to partition; in that case don't break and fall through
                // to SPLIT_EQUAL_COUNTS
                break;
        }
        case SPLIT_EQUAL_COUNTS: {
            // Partition primitives into equally-sized subsets
            mid = (start + end) / 2;
            std::nth_element(&buildData[start], &buildData[mid],
                             &buildData[end-1]+1, ComparePoints(dim));
            break;
        }
        case SPLIT_SAH: default: {
            // Partition primitives using approximate SAH
            if (nPrimitives <= 4) {
                // Partition primitives into equally-sized subsets
                mid = (start + end) / 2;
                std::nth_element(&buildData[start], &buildData[mid],
                                 &buildData[end-1]+1, ComparePoints(dim));
            }
            else {
                // Allocate _BucketInfo_ for SAH partition buckets
                const int nBuckets = 12;
                struct BucketInfo {
                    BucketInfo() { count = 0; }
                    int count;
                    BBox bounds;
                };
                BucketInfo buckets[nBuckets];

                // Initialize _BucketInfo_ for SAH partition buckets
                for (uint32_t i = start; i < end; ++i) {
                    int b = nBuckets *
                        ((buildData[i].centroid[dim] - centroidBounds.pMin[dim]) /
                         (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
                    if (b == nBuckets) b = nBuckets-1;
                    Assert(b >= 0 && b < nBuckets);
                    buckets[b].count++;
                    buckets[b].bounds = Union(buckets[b].bounds, buildData[i].bounds);
                }

                // Compute costs for splitting after each bucket
                float cost[nBuckets-1];
                for (int i = 0; i < nBuckets-1; ++i) {
                    BBox b0, b1;
                    int count0 = 0, count1 = 0;
                    for (int j = 0; j <= i; ++j) {
                        b0 = Union(b0, buckets[j].bounds);
                        count0 += buckets[j].count;
                    }
                    for (int j = i+1; j < nBuckets; ++j) {
                        b1 = Union(b1, buckets[j].bounds);
                        count1 += buckets[j].count;
                    }
                    cost[i] = .125f + (count0*b0.SurfaceArea() + count1*b1.SurfaceArea()) /
                              bbox.SurfaceArea();
                }

                // Find bucket to split at that minimizes SAH metric
                float minCost = cost[0];
                uint32_t minCostSplit = 0;
                for (int i = 1; i < nBuckets-1; ++i) {
                    if (cost[i] < minCost) {
                        minCost = cost[i];
                        minCostSplit = i;
                    }
                }

                // Either create leaf or split primitives at selected SAH bucket
                if (nPrimitives > maxPrimsInNode ||
                    minCost < nPrimitives) {
                    BVHPrimitiveInfo *pmid = std::partition(&buildData[start],
                        &buildData[end-1]+1,
                        CompareToBucket(minCostSplit, nBuckets, dim, centroidBounds));
                    mid = pmid - &buildData[0];
                }
                
                else {
                    // Create leaf _BVHBuildNode_
                    uint32_t firstPrimOffset = orderedPrims.size();
                    for (uint32_t i = start; i < end; ++i) {
                        uint32_t primNum = buildData[i].primitiveNumber;
                        orderedPrims.push_back(primitives[primNum]);
                    }
                    node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
                    return node;
                }
            }
            break;
        }
        }
        node->InitInterior(dim,
                           recursiveBuild(buildArena, buildData, start, mid,
                                          totalNodes, orderedPrims),
                           recursiveBuild(buildArena, buildData, mid, end,
                                          totalNodes, orderedPrims));
    }
    return node;
}


uint32_t BVHAccel::flattenBVHTree(BVHBuildNode *node, uint32_t *offset) {
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    uint32_t myOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        Assert(!node->children[0] && !node->children[1]);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    }
    else {
        // Creater interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVHTree(node->children[0], offset);
        linearNode->secondChildOffset = flattenBVHTree(node->children[1],
                                                       offset);
    }
    return myOffset;
}


BVHAccel::~BVHAccel() {
    FreeAligned(nodes);
}


bool BVHAccel::Intersect(const Ray &ray, Intersection *isect) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTION_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray));
    bool hit = false;
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    // Follow ray through BVH nodes to find primitive intersections
    uint32_t todoOffset = 0, nodeNum = 0;
    uint32_t todo[64];
    while (true) {
        const LinearBVHNode *node = &nodes[nodeNum];
        // Check ray against BVH node
        if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                PBRT_BVH_INTERSECTION_TRAVERSED_LEAF_NODE(const_cast<LinearBVHNode *>(node));
                for (uint32_t i = 0; i < node->nPrimitives; ++i)
                {
                    PBRT_BVH_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                    if (primitives[node->primitivesOffset+i]->Intersect(ray, isect))
                    {
                        PBRT_BVH_INTERSECTION_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                        hit = true;
                    }
                    else {
                        PBRT_BVH_INTERSECTION_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                   }
                }
                if (todoOffset == 0) break;
                nodeNum = todo[--todoOffset];
            }
            else {
                // Put far BVH node on _todo_ stack, advance to near node
                PBRT_BVH_INTERSECTION_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                   todo[todoOffset++] = nodeNum + 1;
                   nodeNum = node->secondChildOffset;
                }
                else {
                   todo[todoOffset++] = node->secondChildOffset;
                   nodeNum = nodeNum + 1;
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            nodeNum = todo[--todoOffset];
        }
    }
    PBRT_BVH_INTERSECTION_FINISHED();
    return hit;
}


bool BVHAccel::IntersectP(const Ray &ray) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTIONP_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray));
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    uint32_t todo[64];
    uint32_t todoOffset = 0, nodeNum = 0;
    while (true) {
        const LinearBVHNode *node = &nodes[nodeNum];
        if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
            // Process BVH node _node_ for traversal
            if (node->nPrimitives > 0) {
                PBRT_BVH_INTERSECTIONP_TRAVERSED_LEAF_NODE(const_cast<LinearBVHNode *>(node));
                  for (uint32_t i = 0; i < node->nPrimitives; ++i) {
                    PBRT_BVH_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
                    if (primitives[node->primitivesOffset+i]->IntersectP(ray)) {
                        PBRT_BVH_INTERSECTIONP_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                        return true;
                    }
                else {
                        PBRT_BVH_INTERSECTIONP_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
                    }
                }
                if (todoOffset == 0) break;
                nodeNum = todo[--todoOffset];
            }
            else {
                PBRT_BVH_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                   /// second child first
                   todo[todoOffset++] = nodeNum + 1;
                   nodeNum = node->secondChildOffset;
                }
                else {
                   todo[todoOffset++] = node->secondChildOffset;
                   nodeNum = nodeNum + 1;
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            nodeNum = todo[--todoOffset];
        }
    }
    PBRT_BVH_INTERSECTIONP_FINISHED();
    return false;
}


BVHAccel *CreateBVHAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps) {
    string splitMethod = ps.FindOneString("splitmethod", "sah");
    uint32_t maxPrimsInNode = ps.FindOneInt("maxnodeprims", 4);
    return new BVHAccel(prims, maxPrimsInNode, splitMethod);
}


