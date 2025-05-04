#include "Scene.h"

Scene::Scene(Camera* p_camera, DXContext* context, CommandListID renderID, DescriptorHeap* dHeap)
	: camera(p_camera),
	pbmpmRP("PBMPMVertexShader.cso", "PBMPMPixelShader.cso", "PBMPMRootSignature.cso", *context, renderID, dHeap),
	pbmpmScene(context, &pbmpmRP, renderToggles),
	//objectRPWire("VertexShader.cso", "PixelShader.cso", "RootSignature.cso", *context, renderID, dHeap),
	//objectSceneGrid(context, &objectRPWire, pbmpmScene.getSimShapes(), 1), 
	//objectSceneSpawners(context, &objectRPWire, pbmpmScene.getSimShapes(), 2), 
	objectRPSolid("VertexShader.cso", "PixelShader.cso", "RootSignature.cso", *context, renderID, dHeap),
	objectSceneSolid(context, &objectRPSolid, pbmpmScene.getSimShapes(), 0),
	// Fluid Mesh Shader Pipeline Construction
	
	fluidRP("VertexShader.cso", "PixelShader.cso", "RootSignature.cso", *context, renderID, dHeap),
	fluidBilevelUniformGridCP("BilevelUniformGridRootSig.cso", "BilevelUniformGrid.cso", *context, CommandListID::FLUID_BILEVEL_UNIFORM_GRID_COMPUTE_ID, dHeap),
	fluidSurfaceBlockDetectionCP("SurfaceBlockDetectionRootSig.cso", "SurfaceBlockDetection.cso", *context, CommandListID::FLUID_SURFACE_BLOCK_DETECTION_COMPUTE_ID, dHeap),
	fluidSurfaceCellDetectionCP("SurfaceCellDetectionRootSig.cso", "SurfaceCellDetection.cso", *context, CommandListID::FLUID_SURFACE_CELL_DETECTION_COMPUTE_ID, dHeap),
	fluidSurfaceVertexCompactionCP("SurfaceVertexCompactionRootSig.cso", "SurfaceVertexCompaction.cso", *context, CommandListID::FLUID_SURFACE_VERTEX_COMPACTION_COMPUTE_ID, dHeap),
	fluidSurfaceVertexDensityCP("SurfaceVertexDensityRootSig.cso", "SurfaceVertexDensity.cso", *context, CommandListID::FLUID_SURFACE_VERTEX_DENSITY_COMPUTE_ID, dHeap),
	fluidSurfaceVertexNormalCP("SurfaceVertexNormalsRootSig.cso", "SurfaceVertexNormals.cso", *context, CommandListID::FLUID_SURFACE_VERTEX_NORMAL_COMPUTE_ID, dHeap),
	fluidMeshPipeline("ConstructMeshShader.cso", "ConstructSurfaceShader.cso", "ConstructMeshRootSig.cso", *context, renderID, dHeap),
	fluidBufferClearCP("bufferClearRootSignature.cso", "bufferClearComputeShader.cso", *context, CommandListID::FLUID_BUFFER_CLEAR_COMPUTE_ID, dHeap),
	fluidDispatchArgDivideCP("DispatchArgDivideRootSig.cso", "DispatchArgDivide.cso", *context, CommandListID::FLUID_DISPATCH_ARG_DIVIDE_COMPUTE_ID, dHeap),
	fluidScene(context, &fluidRP, &fluidBilevelUniformGridCP, &fluidSurfaceBlockDetectionCP, &fluidSurfaceCellDetectionCP, &fluidSurfaceVertexCompactionCP, 
		&fluidSurfaceVertexDensityCP, &fluidSurfaceVertexNormalCP, &fluidBufferClearCP, &fluidDispatchArgDivideCP, &fluidMeshPipeline, 0, 0.488, 1.606, 1.709),

	// Elastic Mesh Shader Pipeline Construction
	elasticRP("VertexShader.cso", "PixelShader.cso", "RootSignature.cso", *context, renderID, dHeap),
	elasticBilevelUniformGridCP("BilevelUniformGridRootSig.cso", "BilevelUniformGrid.cso", *context, CommandListID::ELASTIC_BILEVEL_UNIFORM_GRID_COMPUTE_ID, dHeap),
	elasticSurfaceBlockDetectionCP("SurfaceBlockDetectionRootSig.cso", "SurfaceBlockDetection.cso", *context, CommandListID::ELASTIC_SURFACE_BLOCK_DETECTION_COMPUTE_ID, dHeap),
	elasticSurfaceCellDetectionCP("SurfaceCellDetectionRootSig.cso", "SurfaceCellDetection.cso", *context, CommandListID::ELASTIC_SURFACE_CELL_DETECTION_COMPUTE_ID, dHeap),
	elasticSurfaceVertexCompactionCP("SurfaceVertexCompactionRootSig.cso", "SurfaceVertexCompaction.cso", *context, CommandListID::ELASTIC_SURFACE_VERTEX_COMPACTION_COMPUTE_ID, dHeap),
	elasticSurfaceVertexDensityCP("SurfaceVertexDensityRootSig.cso", "SurfaceVertexDensity.cso", *context, CommandListID::ELASTIC_SURFACE_VERTEX_DENSITY_COMPUTE_ID, dHeap),
	elasticSurfaceVertexNormalCP("SurfaceVertexNormalsRootSig.cso", "SurfaceVertexNormals.cso", *context, CommandListID::ELASTIC_SURFACE_VERTEX_NORMAL_COMPUTE_ID, dHeap),
	elasticMeshPipeline("ConstructMeshShader.cso", "ConstructSurfaceShader.cso", "ConstructMeshRootSig.cso", *context, CommandListID::ELASTIC_MESH_ID, dHeap),
	elasticBufferClearCP("bufferClearRootSignature.cso", "bufferClearComputeShader.cso", *context, CommandListID::ELASTIC_BUFFER_CLEAR_COMPUTE_ID, dHeap),
	elasticDispatchArgDivideCP("DispatchArgDivideRootSig.cso", "DispatchArgDivide.cso", *context, CommandListID::ELASTIC_DISPATCH_ARG_DIVIDE_COMPUTE_ID, dHeap),
	elasticScene(context, &elasticRP, &elasticBilevelUniformGridCP, &elasticSurfaceBlockDetectionCP, &elasticSurfaceCellDetectionCP, &elasticSurfaceVertexCompactionCP, 
		&elasticSurfaceVertexDensityCP, &elasticSurfaceVertexNormalCP, &elasticBufferClearCP, &elasticDispatchArgDivideCP, &elasticMeshPipeline, 1, 0.010, 7.6, 1.010),

	// Sand Mesh Shader Pipeline Construction
	sandRP("VertexShader.cso", "PixelShader.cso", "RootSignature.cso", *context, renderID, dHeap),
	sandBilevelUniformGridCP("BilevelUniformGridRootSig.cso", "BilevelUniformGrid.cso", *context, CommandListID::SAND_BILEVEL_UNIFORM_GRID_COMPUTE_ID, dHeap),
	sandSurfaceBlockDetectionCP("SurfaceBlockDetectionRootSig.cso", "SurfaceBlockDetection.cso", *context, CommandListID::SAND_SURFACE_BLOCK_DETECTION_COMPUTE_ID, dHeap),
	sandSurfaceCellDetectionCP("SurfaceCellDetectionRootSig.cso", "SurfaceCellDetection.cso", *context, CommandListID::SAND_SURFACE_CELL_DETECTION_COMPUTE_ID, dHeap),
	sandSurfaceVertexCompactionCP("SurfaceVertexCompactionRootSig.cso", "SurfaceVertexCompaction.cso", *context, CommandListID::SAND_SURFACE_VERTEX_COMPACTION_COMPUTE_ID, dHeap),
	sandSurfaceVertexDensityCP("SurfaceVertexDensityRootSig.cso", "SurfaceVertexDensity.cso", *context, CommandListID::SAND_SURFACE_VERTEX_DENSITY_COMPUTE_ID, dHeap),
	sandSurfaceVertexNormalCP("SurfaceVertexNormalsRootSig.cso", "SurfaceVertexNormals.cso", *context, CommandListID::SAND_SURFACE_VERTEX_NORMAL_COMPUTE_ID, dHeap),
	sandMeshPipeline("ConstructMeshShader.cso", "ConstructSurfaceShader.cso", "ConstructMeshRootSig.cso", *context, CommandListID::SAND_MESH_ID, dHeap),
	sandBufferClearCP("bufferClearRootSignature.cso", "bufferClearComputeShader.cso", *context, CommandListID::SAND_BUFFER_CLEAR_COMPUTE_ID, dHeap),
	sandDispatchArgDivideCP("DispatchArgDivideRootSig.cso", "DispatchArgDivide.cso", *context, CommandListID::SAND_DISPATCH_ARG_DIVIDE_COMPUTE_ID, dHeap),
	sandScene(context, &sandRP, &sandBilevelUniformGridCP, &sandSurfaceBlockDetectionCP, &sandSurfaceCellDetectionCP, &sandSurfaceVertexCompactionCP,
		&sandSurfaceVertexDensityCP, &sandSurfaceVertexNormalCP, &sandBufferClearCP, &sandDispatchArgDivideCP, &sandMeshPipeline, 2, 0.010, 5.84, 1.180),

	// Visco Mesh Shader Pipeline Construction
	viscoRP("VertexShader.cso", "PixelShader.cso", "RootSignature.cso", *context, renderID, dHeap),
	viscoBilevelUniformGridCP("BilevelUniformGridRootSig.cso", "BilevelUniformGrid.cso", *context, CommandListID::VISCO_BILEVEL_UNIFORM_GRID_COMPUTE_ID, dHeap),
	viscoSurfaceBlockDetectionCP("SurfaceBlockDetectionRootSig.cso", "SurfaceBlockDetection.cso", *context, CommandListID::VISCO_SURFACE_BLOCK_DETECTION_COMPUTE_ID, dHeap),
	viscoSurfaceCellDetectionCP("SurfaceCellDetectionRootSig.cso", "SurfaceCellDetection.cso", *context, CommandListID::VISCO_SURFACE_CELL_DETECTION_COMPUTE_ID, dHeap),
	viscoSurfaceVertexCompactionCP("SurfaceVertexCompactionRootSig.cso", "SurfaceVertexCompaction.cso", *context, CommandListID::VISCO_SURFACE_VERTEX_COMPACTION_COMPUTE_ID, dHeap),
	viscoSurfaceVertexDensityCP("SurfaceVertexDensityRootSig.cso", "SurfaceVertexDensity.cso", *context, CommandListID::VISCO_SURFACE_VERTEX_DENSITY_COMPUTE_ID, dHeap),
	viscoSurfaceVertexNormalCP("SurfaceVertexNormalsRootSig.cso", "SurfaceVertexNormals.cso", *context, CommandListID::VISCO_SURFACE_VERTEX_NORMAL_COMPUTE_ID, dHeap),
	viscoMeshPipeline("ConstructMeshShader.cso", "ConstructSurfaceShader.cso", "ConstructMeshRootSig.cso", *context, CommandListID::VISCO_MESH_ID, dHeap),
	viscoBufferClearCP("bufferClearRootSignature.cso", "bufferClearComputeShader.cso", *context, CommandListID::VISCO_BUFFER_CLEAR_COMPUTE_ID, dHeap),
	viscoDispatchArgDivideCP("DispatchArgDivideRootSig.cso", "DispatchArgDivide.cso", *context, CommandListID::VISCO_DISPATCH_ARG_DIVIDE_COMPUTE_ID, dHeap),
	viscoScene(context, &viscoRP, &viscoBilevelUniformGridCP, &viscoSurfaceBlockDetectionCP, &viscoSurfaceCellDetectionCP, &viscoSurfaceVertexCompactionCP,
		&viscoSurfaceVertexDensityCP, &viscoSurfaceVertexNormalCP, &viscoBufferClearCP, &viscoDispatchArgDivideCP, &viscoMeshPipeline, 3, 0.010, 4.604, 1.010)
		
	// Snow Mesh Shader Pipeline Construction
	/*snowRP("VertexShader.cso", "PixelShader.cso", "RootSignature.cso", *context, CommandListID::ELASTIC_RENDER_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	snowBilevelUniformGridCP("BilevelUniformGridRootSig.cso", "BilevelUniformGrid.cso", *context, CommandListID::SNOW_BILEVEL_UNIFORM_GRID_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 45, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	snowSurfaceBlockDetectionCP("SurfaceBlockDetectionRootSig.cso", "SurfaceBlockDetection.cso", *context, CommandListID::SNOW_SURFACE_BLOCK_DETECTION_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	snowSurfaceCellDetectionCP("SurfaceCellDetectionRootSig.cso", "SurfaceCellDetection.cso", *context, CommandListID::SNOW_SURFACE_CELL_DETECTION_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	snowSurfaceVertexCompactionCP("SurfaceVertexCompactionRootSig.cso", "SurfaceVertexCompaction.cso", *context, CommandListID::SNOW_SURFACE_VERTEX_COMPACTION_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	snowSurfaceVertexDensityCP("SurfaceVertexDensityRootSig.cso", "SurfaceVertexDensity.cso", *context, CommandListID::SNOW_SURFACE_VERTEX_DENSITY_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	snowSurfaceVertexNormalCP("SurfaceVertexNormalsRootSig.cso", "SurfaceVertexNormals.cso", *context, CommandListID::SNOW_SURFACE_VERTEX_NORMAL_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	snowMeshPipeline("ConstructMeshShader.cso", "ConstructSurfaceShader.cso", "ConstructMeshRootSig.cso", *context, CommandListID::SNOW_MESH_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	snowBufferClearCP("bufferClearRootSignature.cso", "bufferClearComputeShader.cso", *context, CommandListID::SNOW_BUFFER_CLEAR_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	snowDispatchArgDivideCP("DispatchArgDivideRootSig.cso", "DispatchArgDivide.cso", *context, CommandListID::SNOW_DISPATCH_ARG_DIVIDE_COMPUTE_ID,
		D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE),
	snowScene(context, &snowRP, &snowBilevelUniformGridCP, &snowSurfaceBlockDetectionCP, &snowSurfaceCellDetectionCP, &snowSurfaceVertexCompactionCP,
		&snowSurfaceVertexDensityCP, &snowSurfaceVertexNormalCP, &snowBufferClearCP, &snowDispatchArgDivideCP, &snowMeshPipeline, 4, 0.010, 7.6, 1.010),*/
{}


//RenderPipeline* Scene::getObjectWirePipeline() {
//	return &objectRPWire;
//}

RenderPipeline* Scene::getObjectSolidPipeline() {
	return &objectRPSolid;
}

RenderPipeline* Scene::getPBMPMRenderPipeline() {
	return &pbmpmRP;
}

MeshPipeline* Scene::getFluidMeshPipeline() {
	return &fluidMeshPipeline;
}

MeshPipeline* Scene::getElasticMeshPipeline() {
	return &elasticMeshPipeline;
}

MeshPipeline* Scene::getSandMeshPipeline() {
	return &sandMeshPipeline;
}

MeshPipeline* Scene::getViscoMeshPipeline() {
	return &viscoMeshPipeline;
}

//MeshPipeline* Scene::getSnowMeshPipeline() {
//	return &snowMeshPipeline;
//}
	
void Scene::compute(float isMeshShading) {
	pbmpmScene.compute();
	int particles = pbmpmScene.transferAndGetNumParticles();
	if (isMeshShading) {
		if (renderToggles[0]) {
			fluidScene.compute(
				pbmpmScene.getPositionBuffer(),
				particles
			);
		}
		if (renderToggles[1]) {
			elasticScene.compute(
				pbmpmScene.getPositionBuffer(),
				particles
			);
		}
		if (renderToggles[2]) {
			sandScene.compute(
				pbmpmScene.getPositionBuffer(),
				particles
			);
		}
		if (renderToggles[3]) {
			viscoScene.compute(
				pbmpmScene.getPositionBuffer(),
				particles
			);
		}
		//if (renderToggles[4]) {
		//	snowScene.compute(
		//		pbmpmScene.getPositionBuffer(),
		//		particles
		//	);
		//}
	}
}

void Scene::drawPBMPM() {
	pbmpmScene.draw(camera);
}

void Scene::drawFluid(unsigned int renderMeshlets, unsigned int renderOptions) {
	fluidScene.draw(camera, renderMeshlets, renderOptions);
}

void Scene::drawElastic(unsigned int renderMeshlets, unsigned int renderOptions) {
	elasticScene.draw(camera, renderMeshlets, renderOptions);
}

void Scene::drawSand(unsigned int renderMeshlets, unsigned int renderOptions) {
	sandScene.draw(camera, renderMeshlets, renderOptions);
}

void Scene::drawVisco(unsigned int renderMeshlets, unsigned int renderOptions) {
	viscoScene.draw(camera, renderMeshlets, renderOptions);
}

//void Scene::drawSnow(unsigned int renderMeshlets, unsigned int renderOptions) {
//	snowScene.draw(camera, renderMeshlets, renderOptions);
//}

//void Scene::drawGrid() {
//	objectSceneGrid.draw(camera);
//}
//
//void Scene::drawSpawners() {
//	objectSceneSpawners.draw(camera);
//}

void Scene::drawSolidObjects(XMFLOAT3& leftPos, XMVECTOR& leftRot, XMFLOAT3& rightPos, XMVECTOR& rightRot) {
	objectSceneSolid.draw(camera, leftPos, leftRot, rightPos, rightRot);
}

void Scene::releaseResources() {
	//objectSceneGrid.releaseResources();
	//objectSceneSpawners.releaseResources();
	objectSceneSolid.releaseResources();
	pbmpmScene.releaseResources();
	fluidScene.releaseResources();
	elasticScene.releaseResources();
	viscoScene.releaseResources();
	sandScene.releaseResources();
	pbmpmRP.releaseResources();
	fluidRP.releaseResources();
	elasticRP.releaseResources();
	sandRP.releaseResources();
	viscoRP.releaseResources();
	//snowRP.releaseResources();
	//snowScene.releaseResources();
}