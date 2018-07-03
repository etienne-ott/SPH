#include "domain.h"
#include <stdexcept>
#include <string>
#include <cmath>
#include "util/random_pool.h"

Domain::Domain(DomainType type, YAML::Node& params) {
	this->dtype = type;
	this->size = new float[3];
	this->offset = new float[3];
	this->boundingBox = new float[6];

	float tmp[3] = {
		params["offset_x"].as<float>(),
		params["offset_y"].as<float>(),
		params["offset_z"].as<float>()
	};
	this->setOffset(tmp);

	float tmp2[3] = {
		params["size_x"].as<float>(),
		params["size_y"].as<float>(),
		params["size_z"].as<float>()
	};
	this->setSize(tmp2);
}

Domain::~Domain() {
	delete[] size;
	delete[] offset;
	delete[] boundingBox;
}

void Domain::setSize(float* newSize) {
	for (int i = 0; i < 3; i++) {
		this->size[i] = newSize[i];
	}

	this->setBoundingBox();
}

void Domain::setOffset(float* off) {
	for (int i = 0; i < 3; i++) {
		this->offset[i] = off[i];
	}

	this->setBoundingBox();
}

void Domain::setBoundingBox() {
	if (this->dtype == DomainType::Cube) {
		for (int i = 0; i < 3; i++) {
			this->boundingBox[i] = this->offset[i];
			this->boundingBox[i + 3] = this->offset[i] + this->size[i];
		}
	} else if (this->dtype == DomainType::Sphere) {
		for (int i = 0; i < 3; i++) {
			this->boundingBox[i] = this->offset[i] - this->size[i];
			this->boundingBox[i + 3] = this->offset[i] + this->size[i];
		}
	} else if (this->dtype == DomainType::Mesh) {
		// Nothing to do, this is handled in setMesh
	} else {
		throw std::runtime_error(std::string("Trying to establish boundary box of a domain with unspecified type."));
	}
}

float* Domain::getBoundingBox() {
	return this->boundingBox;
}

float* Domain::getSize() {
	return this->size;
}

float* Domain::getOffset() {
	return this->offset;
}

DomainType Domain::getType() {
	return this->dtype;
}

void Domain::setMesh(Mesh* mesh) {
	this->mesh = mesh;
	float* box = this->mesh->getBoundingBox();
	this->boundingBox[0] = box[0];
	this->boundingBox[1] = box[1];
	this->boundingBox[2] = box[2];
	this->boundingBox[3] = box[3];
	this->boundingBox[4] = box[4];
	this->boundingBox[5] = box[5];
}

Mesh* Domain::getMesh() {
	return this->mesh;
}

float Domain::getVolume() {
	if (this->dtype == DomainType::Cube) {
		return this->size[0] * this->size[1] * this->size[2];

	} else if (this->dtype == DomainType::Sphere) {
		return 4.0 * M_PI * this->size[0] * this->size[1] * this->size[2] / 3.0;

	} else if (this->dtype == DomainType::Mesh) {
		return this->mesh->getVolume();

	} else {
		throw std::runtime_error(std::string("Calling getVolume on a domain with undefined type."));
	}
}

bool Domain::pointIsInDomain(float* position) {
	if (this->dtype == DomainType::Cube) {
		bool pointWithin = true;

		for (int i = 0; i < 3; i++) {
			pointWithin = pointWithin
				&& position[i] >= this->offset[i]
				&& position[i] <= this->offset[i] + this->size[i];
		}

		return pointWithin;

	} else if (this->dtype == DomainType::Sphere) {
		float sum = 0.0;
		float tmp;
		for (int i = 0; i < 3; i++) {
			tmp = (position[i] - this->offset[i]) / this->size[i];
			sum += tmp * tmp;
		}

		return sum <= 1.0;

	} else if (this->dtype == DomainType::Mesh) {
		Vector3D<float> pos = Vector3D<float>(position[0], position[1], position[2]);
		return this->mesh->pointIsInsideMesh(pos);

	} else {
		throw std::runtime_error(std::string("Calling pointIsInDomain on a domain with undefined type."));
	}
}
