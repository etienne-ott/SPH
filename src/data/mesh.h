#pragma once
#include <vector>
#include <string>
#include "util/random_pool.h"
#include "data/vector3D.h"

#define NR_RECALC_FLAGS 3
enum RecalculationFlags : int {
    BoundingBox = 0,
    FaceNormals = 1,
    Volume = 2
};

class Mesh
{
protected:
    std::vector<int> faces;
    std::vector<Vector3D<float>> vertices;
    std::vector<Vector3D<float>> faceNormals;
    float* boundingBox;
    float volume;
    bool* needsRecalculation;

    void calculateFaceNormals();
    void calculateBoundingBox();
    void calculateVolume();

public:
    Mesh();
    ~Mesh();

    std::vector<int>& getFaces() {return this->faces;}
    std::vector<Vector3D<float>>& getVertices() {return this->vertices;}

    float getVolume();
    float* getBoundingBox();
    std::vector<Vector3D<float>>& getFaceNormals();

    void loadMeshFromOBJFile(std::string filepath);
    void writeMeshToOBJFile(std::string filepath);

    void centerOnOrigin();
    void centerOn(Vector3D<float> point);
    void scaleToNormal();
    void scaleTo(float scale);
    void rotate(float* matrix);

    bool pointIsInsideMesh(Vector3D<float>& x);
    bool rayIntersectsFace(Vector3D<float>& x, Vector3D<float>& r, int faceIdx);
    bool lineIntersectsFace(
        Vector3D<float>& start,
        Vector3D<float>& line,
        int faceIdx,
        Vector3D<float>& intersection
    );
};