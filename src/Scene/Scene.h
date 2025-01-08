#pragma once

#include "ObjectScene.h"
#include "PBMPMScene.h"
#include "MeshShadingScene.h"
#include "../D3D/Pipeline/RenderPipeline.h"
#include "../D3D/Pipeline/ComputePipeline.h"
#include "../D3D/Pipeline/MeshPipeline.h"

class Scene {
public:
	Scene() = delete;
	Scene(Camera* camera, DXContext* context);

	RenderPipeline* getObjectWirePipeline();
	RenderPipeline* getObjectSolidPipeline();
	RenderPipeline* getPBMPMRenderPipeline();
	MeshPipeline* getFluidMeshPipeline();

	void compute(float isMeshShading = false);
	void drawPBMPM();
	void drawFluid(unsigned int renderMeshlets, unsigned int renderOptions);
	void drawGrid();
	void drawSpawners();
	void drawSolidObjects();

	void releaseResources();

	PBMPMConstants getPBMPMConstants() { return pbmpmScene.getConstants(); }
	void updatePBMPMConstants(PBMPMConstants& newConstants);

	float* getFluidIsovalue() { return fluidScene.getIsovalue(); }
	float* getFluidKernelScale() { return fluidScene.getKernelScale(); }
	float* getFluidKernelRadius() { return fluidScene.getKernelRadius(); }

	unsigned int* getPBMPMSubstepCount() { return pbmpmScene.getSubstepCount(); }

	int getNumParticles() { return pbmpmScene.getNumParticles(); }

private:
	Camera* camera;

	RenderPipeline pbmpmRP;
	PBMPMScene pbmpmScene;

	RenderPipeline objectRPWire;
	ObjectScene objectSceneGrid;
	ObjectScene objectSceneSpawners;
	RenderPipeline objectRPSolid;
	ObjectScene objectSceneSolid;

	// Fluid Mesh
	RenderPipeline fluidRP;
	ComputePipeline fluidBilevelUniformGridCP;
	ComputePipeline fluidSurfaceBlockDetectionCP;
	ComputePipeline fluidSurfaceCellDetectionCP;
	ComputePipeline fluidSurfaceVertexCompactionCP;
	ComputePipeline fluidSurfaceVertexDensityCP;
	ComputePipeline fluidSurfaceVertexNormalCP;
	ComputePipeline fluidBufferClearCP;
	MeshPipeline fluidMeshPipeline;
	MeshShadingScene fluidScene;

	RenderPipeline* currentRP;
	ComputePipeline* currentCP;
};
