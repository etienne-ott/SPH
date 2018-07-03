#pragma once

enum class DomainType : int {
	Cube = 1,
	Sphere = 2,
	Mesh = 3
};

enum class DistributionType : int {
	VolumeGrid = 1,
	WhiteNoise = 2,
	FastPoissonDisk = 3,
	GoldenSet = 4,
	Hammersley = 5,
	Halton = 6,
	SpherePacking = 7
};