#include "data/mesh.h"
#include <cstring>
#include <cmath>
#include "util/misc_math.h"
#include "util/random_pool.h"

Mesh::Mesh()
{
    this->boundingBox = new float[6];
    this->faces = std::vector<int>();
    this->vertices = std::vector<Vector3D<float>>();
    this->faceNormals = std::vector<Vector3D<float>>();
    this->needsRecalculation = new bool[NR_RECALC_FLAGS];

    for (uint i = 0; i < NR_RECALC_FLAGS; i++) {
        this->needsRecalculation[i] = true;
    }
}

Mesh::~Mesh()
{
    delete[] this->boundingBox;
    delete[] this->needsRecalculation;
}

void Mesh::calculateFaceNormals()
{
    this->faceNormals.clear();

    uint k = 0;

    for (uint i = 0; i * 3 < this->faces.size(); i++) {
        k = 3 * i;

        Vector3D<float> a = this->vertices[this->faces[k]];
        Vector3D<float> b = this->vertices[this->faces[k + 1]];
        Vector3D<float> c = this->vertices[this->faces[k + 2]];

        Vector3D<float> tmp1 = (b - a) % (c - a);
        Vector3D<float>* tmp2 = new Vector3D<float>(tmp1 / tmp1.getLen());

        this->faceNormals.push_back(*tmp2);
    }

    this->needsRecalculation[RecalculationFlags::FaceNormals] = false;
}

void Mesh::loadMeshFromOBJFile(std::string filepath)
{
    FILE* file = fopen(filepath.c_str(), "r");
    if (file == NULL) {
        printf("Can't open mesh file %s!\n", filepath.c_str());
        return;
    }

    bool foundZeroIndex = false;
    int result;

    while (true) {
        char lineBuffer[1024];

        result = fscanf(file, "%s", lineBuffer);
        if (result == EOF) {
            break;
        }

        if (std::strcmp(lineBuffer, "v") == 0) {
            float coords[3];
            result = fscanf(file, "%f %f %f\n", &coords[0], &coords[1], &coords[2]);
            this->vertices.push_back(Vector3D<float>(coords));

        } else if (std::strcmp(lineBuffer, "vt") == 0) {
            // we ignore textures
            result = fscanf(file, "%s\n", lineBuffer);

        } else if (std::strcmp(lineBuffer, "vn") == 0 ) {
            // we ignore vertex normals
            result = fscanf(file, "%s\n", lineBuffer);

        } else if (std::strcmp(lineBuffer, "f") == 0 ) {
            int vertexIndex[3], uvIndex[3], normalIndex[3];
            bool foundMatch = false;

            // load line into char buffer so we can read it multiple times
            // this syntax is not a regex pattern, it's scanf's weird formating
            // (although it is technically a regular expression)
            result = fscanf(file, "%[^\n]s", lineBuffer);

            // scan for triple vertex/texture/normal format
            int matches = sscanf(
                lineBuffer, "%d/%d/%d %d/%d/%d %d/%d/%d",
                &vertexIndex[0], &uvIndex[0], &normalIndex[0],
                &vertexIndex[1], &uvIndex[1], &normalIndex[1],
                &vertexIndex[2], &uvIndex[2], &normalIndex[2]
            );

            if (matches == 9) {
                foundMatch = true;
            }

            // scan for triple vertex/texture format
            if (!foundMatch) {
                matches = sscanf(
                    lineBuffer, "%d/%d %d/%d %d/%d",
                    &vertexIndex[0], &uvIndex[0],
                    &vertexIndex[1], &uvIndex[1],
                    &vertexIndex[2], &uvIndex[2]
                );

                if (matches == 6) {
                    foundMatch = true;
                }
            }

            // scan for triple vertex format
            if (!foundMatch) {
                matches = sscanf(
                    lineBuffer, "%d %d %d",
                    &vertexIndex[0],
                    &vertexIndex[1],
                    &vertexIndex[2]
                );

                if (matches == 3) {
                    foundMatch = true;
                }
            }

            if (!foundMatch) {
                // Okay, we have no idea what format this is
                printf("Mesh file format not supported. More info in documentation.\n");
                fclose(file);
                return;
            }

            // @TODO handle textures and vertex normals too
            this->faces.push_back(vertexIndex[0]);
            this->faces.push_back(vertexIndex[1]);
            this->faces.push_back(vertexIndex[2]);

            // check if we need to correct indizes down by one later on
            if (vertexIndex[0] == 0 || vertexIndex[1] == 0 || vertexIndex[2] == 0) {
                foundZeroIndex = true;
            }

        } else {
            // other line, ignore it entirely
            result = fscanf(file, "%s\n", lineBuffer);
        }
    }

    fclose(file);

    // correct indizes if necessary
    if (!foundZeroIndex) {
        for (uint i = 0; i < this->faces.size(); i++) {
            this->faces[i] -= 1;
        }
    }

    this->needsRecalculation[RecalculationFlags::BoundingBox] = true;
    this->needsRecalculation[RecalculationFlags::FaceNormals] = true;
    this->needsRecalculation[RecalculationFlags::Volume] = true;
}

void Mesh::writeMeshToOBJFile(std::string filepath)
{
    FILE* file = fopen(filepath.c_str(), "w");
    if (file == NULL) {
        printf("Can't open mesh file %s to write!\n", filepath.c_str());
        return;
    }

    for (uint i = 0; i < this->vertices.size(); i++) {
        fprintf(
            file,
            "v %f %f %f\n",
            this->vertices[i].getX(),
            this->vertices[i].getY(),
            this->vertices[i].getZ()
        );
    }

    for (uint i = 0; i < this->faces.size(); i+=3) {
        fprintf(
            file,
            "f %d %d %d\n",
            this->faces[i] + 1,     // +1 because we shifted the indizes
            this->faces[i + 1] + 1, // during loading
            this->faces[i + 2] + 1
        );
    }

    fclose(file);
}

void Mesh::centerOnOrigin()
{
    Vector3D<float> og = Vector3D<float>(0.f, 0.f, 0.f);
    this->centerOn(og);
}

void Mesh::centerOn(Vector3D<float> point) {
    Vector3D<float> cog = Vector3D<float>(0.f, 0.f, 0.f);
    Vector3D<float> centeredVertex;

    for (uint i = 0; i < this->vertices.size(); i++) {
        cog += this->vertices[i];
    }
    cog /= float(this->vertices.size());

    for (uint i = 0; i < this->vertices.size(); i++) {
        centeredVertex = Vector3D<float>(this->vertices[i] - cog + point);
        this->vertices[i] = centeredVertex;
    }

    this->needsRecalculation[RecalculationFlags::BoundingBox] = true;
}

void Mesh::scaleTo(float scale)
{
    float* box = this->getBoundingBox();
    float boxScale = scale / std::max(
        fabs(box[3] - box[0]),
        std::max(
            fabs(box[3] - box[1]),
            fabs(box[5] - box[2])
        )
    );

    for (uint i = 0; i < this->vertices.size(); i++) {
        this->vertices[i].set(
            this->vertices[i].getX() * boxScale,
            this->vertices[i].getY() * boxScale,
            this->vertices[i].getZ() * boxScale
        );
    }

    this->needsRecalculation[RecalculationFlags::BoundingBox] = true;
    this->needsRecalculation[RecalculationFlags::Volume] = true;
}

void Mesh::rotate(float* matrix)
{
    float newCoords[3];

    for (uint i = 0; i < this->vertices.size(); i++) {
        newCoords[0] =
            this->vertices[i].getX() * matrix[0]
            + this->vertices[i].getY() * matrix[1]
            + this->vertices[i].getZ() * matrix[2];
        newCoords[1] =
            this->vertices[i].getX() * matrix[3]
            + this->vertices[i].getY() * matrix[4]
            + this->vertices[i].getZ() * matrix[5];
        newCoords[2] =
            this->vertices[i].getX() * matrix[6]
            + this->vertices[i].getY() * matrix[7]
            + this->vertices[i].getZ() * matrix[8];

        this->vertices[i].setv(newCoords);
    }

    for (uint i = 0; i < this->faceNormals.size(); i++) {
        newCoords[0] =
            this->faceNormals[i].getX() * matrix[0]
            + this->faceNormals[i].getY() * matrix[1]
            + this->faceNormals[i].getZ() * matrix[2];
        newCoords[1] =
            this->faceNormals[i].getX() * matrix[3]
            + this->faceNormals[i].getY() * matrix[4]
            + this->faceNormals[i].getZ() * matrix[5];
        newCoords[2] =
            this->faceNormals[i].getX() * matrix[6]
            + this->faceNormals[i].getY() * matrix[7]
            + this->faceNormals[i].getZ() * matrix[8];

        this->faceNormals[i].setv(newCoords);
    }

    this->needsRecalculation[RecalculationFlags::BoundingBox] = true;
}

bool Mesh::pointIsInsideMesh(Vector3D<float>& x)
{
    // @TODO check if point is on faces first, can be omitted
    // for increased performance

    // now check ray in random direction with length of bounding
    // box diagonal. if no or even number of intersections, the
    // point is on the outside of the mesh
    int intersections = 0;

    // @TODO check if hardcoded seed is fine
    RandomPool pool = RandomPool(long(12345));
    float* box = this->getBoundingBox();
    float rayLength = fastSqrt2(
        (box[3] - box[0]) * (box[3] - box[0])
        + (box[4] - box[1]) * (box[4] - box[1])
        + (box[5] - box[2]) * (box[5] - box[2])
    );

    // @TODO this random direction is not isotropic on a sphere,
    // but isotropic on a cube. Does this matter?
    Vector3D<float> dx = Vector3D<float>(
        pool.nextFloat(0.0, 2.0) * rayLength,
        pool.nextFloat(0.0, 2.0) * rayLength,
        pool.nextFloat(0.0, 2.0) * rayLength
    );

    for (uint i = 0, k = 0; k < this->faces.size(); k = 3 * ++i) {
        if (this->rayIntersectsFace(x, dx, i)) {
            intersections++;
        }
    }

    return intersections > 0 && intersections % 2 == 1;
}

// This is the Möller–Trumbore intersection algorithm
// see https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
bool Mesh::rayIntersectsFace(Vector3D<float>& x, Vector3D<float>& r, int faceIdx)
{
    const float EPSILON = 0.0000001;
    Vector3D<float> vertex0 = this->vertices[this->faces[faceIdx * 3]];
    Vector3D<float> vertex1 = this->vertices[this->faces[faceIdx * 3 + 1]];
    Vector3D<float> vertex2 = this->vertices[this->faces[faceIdx * 3 + 2]];
    Vector3D<float> edge1, edge2, h, s, q;
    float a,f,u,v;

    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;

    h = r % edge2;
    a = h * edge1;

    if (a > -EPSILON && a < EPSILON) {
        return false;
    }

    f = 1.f / a;
    s = x - vertex0;
    u = f * (s * h);

    if (u < 0.0 || u > 1.0) {
        return false;
    }

    q = s % edge1;
    v = f * (r * q);

    if (v < 0.0 || u + v > 1.0) {
        return false;
    }

    // At this stage we can compute t to find out where the intersection point is on the line.
    float t = f * (edge2 * q);

    if (t > EPSILON) {
        // ray intersection
        return true;
    } else {
        // This means that there is a line intersection but not a ray intersection.
        return false;
    }
}

bool Mesh::lineIntersectsFace(
    Vector3D<float>& start,
    Vector3D<float>& line,
    int faceIdx,
    Vector3D<float>& intersection
) {
    const float EPSILON = 0.0000001;
    Vector3D<float> vertex0 = this->vertices[this->faces[faceIdx * 3]];
    Vector3D<float> vertex1 = this->vertices[this->faces[faceIdx * 3 + 1]];
    Vector3D<float> vertex2 = this->vertices[this->faces[faceIdx * 3 + 2]];
    Vector3D<float> edge1, edge2, h, s, q;
    float a,f,u,v;

    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;

    h = line % edge2;
    a = h * edge1;

    if (a > -EPSILON && a < EPSILON) {
        return false;
    }

    f = 1.f / a;
    s = start - vertex0;
    u = f * (s * h);

    if (u < 0.0 || u > 1.0) {
        return false;
    }

    q = s % edge1;
    v = f * (line * q);

    if (v < 0.0 || u + v > 1.0) {
        return false;
    }

    // At this stage we can compute t to find out where the intersection point is on the line.
    float t = f * (edge2 * q);

    if (t > EPSILON) {
        // ray intersection, now check if the intersection is on the line
        // or on the extension (ray) of it
        Vector3D<float> inter = Vector3D<float>(start + line * t);
        intersection.set(inter.getX(), inter.getY(), inter.getZ());
        return t <= 1.f;
    } else {
        // This means that there is a line intersection but not a ray intersection,
        // meaning on the other side of the line
        return false;
    }
}

void Mesh::calculateVolume()
{
    // Volume for a mesh is calculated as the sum of signed volumes
    // of each face, which itself is the dot product of the first
    // vertex with the crossproduct of the other two, divided by 6
    float sum = 0.f;
    for (uint i = 0; i < this->faces.size(); i += 3) {
        sum += this->vertices[this->faces[i]]
        * (
           this->vertices[this->faces[i + 1]]
           % this->vertices[this->faces[i + 2]]
        );
    }

    this->volume = fabs(sum / 6.f);
    this->needsRecalculation[RecalculationFlags::Volume] = false;
}

float Mesh::getVolume()
{
    if (this->needsRecalculation[RecalculationFlags::Volume]) {
        this->calculateVolume();
    }
    return this->volume;
}

std::vector<Vector3D<float>>& Mesh::getFaceNormals()
{
    if (this->needsRecalculation[RecalculationFlags::FaceNormals]) {
        this->calculateFaceNormals();
    }
    return this->faceNormals;
}

float* Mesh::getBoundingBox()
{
    if (this->needsRecalculation[RecalculationFlags::BoundingBox]) {
        this->calculateBoundingBox();
    }
    return this->boundingBox;
}

void Mesh::calculateBoundingBox()
{
    // @TODO something with float.maxVal or something?
    this->boundingBox[0] = 1e9;
    this->boundingBox[1] = 1e9;
    this->boundingBox[2] = 1e9;
    this->boundingBox[3] = -1e9;
    this->boundingBox[4] = -1e9;
    this->boundingBox[5] = -1e9;

    for (uint i = 0; i < this->vertices.size(); i++) {
        if (this->vertices[i].getX() < this->boundingBox[0])
            this->boundingBox[0] = this->vertices[i].getX();
        if (this->vertices[i].getY() < this->boundingBox[1])
            this->boundingBox[1] = this->vertices[i].getY();
        if (this->vertices[i].getZ() < this->boundingBox[2])
            this->boundingBox[2] = this->vertices[i].getZ();
        if (this->vertices[i].getX() > this->boundingBox[3])
            this->boundingBox[3] = this->vertices[i].getX();
        if (this->vertices[i].getY() > this->boundingBox[4])
            this->boundingBox[4] = this->vertices[i].getY();
        if (this->vertices[i].getZ() > this->boundingBox[5])
            this->boundingBox[5] = this->vertices[i].getZ();
    }

    this->needsRecalculation[RecalculationFlags::BoundingBox] = false;
}
