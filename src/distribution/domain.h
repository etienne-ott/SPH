#pragma once

#include "distribution/distributionEnums.h"
#include "data/mesh.h"
#include <yaml-cpp/yaml.h>

class Domain {
private:
	DomainType dtype;
	float* size;
	float* offset;
	float* boundingBox;
	Mesh* mesh;

	void setBoundingBox();

public:
	Domain(DomainType type, YAML::Node& params);
	~Domain();

	void setSize(float* size);
	float* getSize();
	void setOffset(float* off);
	float* getOffset();
	DomainType getType();
	void setMesh(Mesh* mesh);
	Mesh* getMesh();

	float getVolume();

	float* getBoundingBox();

	bool pointIsInDomain(float* position);
};