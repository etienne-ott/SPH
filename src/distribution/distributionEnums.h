#pragma once

enum class DomainType : int {
	Cube = 1,
	Sphere = 2,
	Mesh = 3
};

enum class DistributionType : int {
	VolumeGrid = 1,
	SpherePacking = 2,
	WhiteNoise = 3,
	FastPoissonDisk = 4,
	Hammersley = 5,
	Halton = 6,
	GoldenSet = 7
};