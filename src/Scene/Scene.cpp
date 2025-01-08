#include "Scene.h"

Scene::Scene(Camera* p_camera, DXContext* context)
	:  camera(p_camera),
	pbmpmRP("PBMPMVertexShader.cso", "PBMPMPixelShader.cso", "PBMPMRootSignature.cso", *context, CommandListID::PBMPM_RENDER_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	pbmpmScene(context, &pbmpmRP),
	objectRPWire("VertexShader.cso", "PixelShader.cso", "RootSignature.cso", *context, CommandListID::OBJECT_RENDER_WIRE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	objectSceneGrid(context, &objectRPWire, pbmpmScene.getSimShapes(), 1), 
	objectSceneSpawners(context, &objectRPWire, pbmpmScene.getSimShapes(), 2), 
	objectRPSolid("VertexShader.cso", "PixelShader.cso", "RootSignature.cso", *context, CommandListID::OBJECT_RENDER_SOLID_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	objectSceneSolid(context, &objectRPSolid, pbmpmScene.getSimShapes(), 0),
	fluidRP("VertexShader.cso", "PixelShader.cso", "RootSignature.cso", *context, CommandListID::FLUID_RENDER_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	fluidBilevelUniformGridCP("BilevelUniformGridRootSig.cso", "BilevelUniformGrid.cso", *context, CommandListID::BILEVEL_UNIFORM_GRID_COMPUTE_ID, 
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 45, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	fluidSurfaceBlockDetectionCP("SurfaceBlockDetectionRootSig.cso", "SurfaceBlockDetection.cso", *context, CommandListID::SURFACE_BLOCK_DETECTION_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	fluidSurfaceCellDetectionCP("SurfaceCellDetectionRootSig.cso", "SurfaceCellDetection.cso", *context, CommandListID::SURFACE_CELL_DETECTION_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	fluidSurfaceVertexCompactionCP("SurfaceVertexCompactionRootSig.cso", "SurfaceVertexCompaction.cso", *context, CommandListID::SURFACE_VERTEX_COMPACTION_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	fluidSurfaceVertexDensityCP("SurfaceVertexDensityRootSig.cso", "SurfaceVertexDensity.cso", *context, CommandListID::SURFACE_VERTEX_DENSITY_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	fluidSurfaceVertexNormalCP("SurfaceVertexNormalsRootSig.cso", "SurfaceVertexNormals.cso", *context, CommandListID::SURFACE_VERTEX_NORMAL_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	fluidMeshPipeline("ConstructMeshShader.cso", "ConstructSurfaceShader.cso", "ConstructMeshRootSig.cso", *context, CommandListID::FLUID_MESH_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	fluidBufferClearCP("bufferClearRootSignature.cso", "bufferClearComputeShader.cso", *context, CommandListID::FLUID_BUFFER_CLEAR_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	fluidScene(context, &fluidRP, &fluidBilevelUniformGridCP, &fluidSurfaceBlockDetectionCP, &fluidSurfaceCellDetectionCP, &fluidSurfaceVertexCompactionCP, &fluidSurfaceVertexDensityCP, &fluidSurfaceVertexNormalCP, &fluidBufferClearCP, &fluidMeshPipeline),
	currentRP(),
	currentCP()
{}

RenderPipeline* Scene::getObjectWirePipeline() {
	return &objectRPWire;
}

RenderPipeline* Scene::getObjectSolidPipeline() {
	return &objectRPSolid;
}

RenderPipeline* Scene::getPBMPMRenderPipeline() {
	return &pbmpmRP;
}

MeshPipeline* Scene::getFluidMeshPipeline() {
	return &fluidMeshPipeline;
}

void Scene::compute(float isMeshShading) {
	pbmpmScene.compute();
	int particles = pbmpmScene.transferAndGetNumParticles();
	if (isMeshShading) {
		fluidScene.compute(
			pbmpmScene.getPositionBuffer(),
			particles
		);
	}
}

void Scene::drawPBMPM() {
	pbmpmScene.draw(camera);
}

void Scene::drawFluid(unsigned int renderMeshlets, unsigned int renderOptions) {
	fluidScene.draw(camera, renderMeshlets, renderOptions);
}

void Scene::drawGrid() {
	objectSceneGrid.draw(camera);
}

void Scene::drawSpawners() {
	objectSceneSpawners.draw(camera);
}

void Scene::drawSolidObjects() {
	objectSceneSolid.draw(camera);
}

void Scene::releaseResources() {
	objectSceneGrid.releaseResources();
	objectSceneSpawners.releaseResources();
	objectSceneSolid.releaseResources();
	pbmpmScene.releaseResources();
	fluidScene.releaseResources();
}

void Scene::updatePBMPMConstants(PBMPMConstants& newConstants) {
	pbmpmScene.updateConstants(newConstants);
}