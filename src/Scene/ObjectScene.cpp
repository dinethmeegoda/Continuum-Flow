#include "ObjectScene.h"
#include "SceneConstants.h"

ObjectScene::ObjectScene(DXContext* context, RenderPipeline* pipeline, std::vector<SimShape>& shapes, int renderWireframe)
	: Drawable(context, pipeline), shapes(shapes)
{
    if (renderWireframe == 1) {
        constructSceneGrid();
    }
	else if (renderWireframe == 2) {
		constructSceneSpawners();
	}
    else {
        constructSceneSolid();
    }
}

void ObjectScene::constructSceneGrid() {
	renderPipeline->createPSOD();
	renderPipeline->createPipelineState(context->getDevice());

    std::vector<std::string> inputStrings;
    inputStrings.push_back("objs\\cube.obj");

    //cube for grid
    XMFLOAT4X4 gridModelMatrix;
    XMStoreFloat4x4(&gridModelMatrix, XMMatrixScaling(GRID_WIDTH, GRID_HEIGHT, GRID_DEPTH));
    modelMatrices.push_back(gridModelMatrix);

	// vector for colors of grid lines
	XMFLOAT3 color = XMFLOAT3(0.0f, 1.0f, 0.0f);

    auto string = inputStrings.front();
    auto m = modelMatrices.front();
    Mesh newMesh = Mesh((std::filesystem::current_path() / string).string(), context, renderPipeline->getCommandList(), renderPipeline, m, true, color);
    meshes.push_back(newMesh);
    sceneSize += newMesh.getNumTriangles();
}

void ObjectScene::constructSceneSpawners() {
    renderPipeline->createPSOD();
    renderPipeline->createPipelineState(context->getDevice());

    std::vector<std::string> inputStrings;

    // vector for colors of lines
    std::vector<XMFLOAT3> colors = {};

    //cubes for shapes
    for (const auto& shape : shapes) {
        inputStrings.push_back("objs\\cube.obj");

        // Decide color based on shape function
        XMFLOAT3 color;
        switch (shape.functionality) {
        case 0: // Emitter
            color = XMFLOAT3(0.0f, 0.7f, 0.8f);
            break;
        case 1: // Obstacle
            color = XMFLOAT3(0.0f, 0.0f, 0.0f);
            break;
        case 2: // Drain
            color = XMFLOAT3(1.0f, 0.0f, 0.25f);
            break;
        case 3: // Initial Emitter
            color = XMFLOAT3(0.0f, 0.0f, 0.9f);
            break;
        default:
            color = XMFLOAT3(1.0f, 1.0f, 1.0f);
            break;
        }
        colors.push_back(color);

        XMFLOAT4X4 simShapeMatrix;
        XMStoreFloat4x4(&simShapeMatrix, XMMatrixMultiply(
            XMMatrixScaling(shape.halfSize.x * 2, shape.halfSize.y * 2, shape.halfSize.z * 2),
            XMMatrixTranslation(shape.position.x - shape.halfSize.x, shape.position.y - shape.halfSize.y, shape.position.z - shape.halfSize.z)
        ));
        modelMatrices.push_back(simShapeMatrix);
    }

    //push shapes as wireframe
    for (int i = 0; i < inputStrings.size(); i++) {
        auto string = inputStrings.at(i);
        auto m = modelMatrices.at(i);
        Mesh newMesh = Mesh((std::filesystem::current_path() / string).string(), context, renderPipeline->getCommandList(), renderPipeline, m, true, colors.at(i));
        meshes.push_back(newMesh);
        sceneSize += newMesh.getNumTriangles();
    }
}

void ObjectScene::constructSceneSolid() {
    //cube for ground
    std::vector<std::string> inputStrings;
    inputStrings.push_back("objs\\cube.obj");

    XMFLOAT4X4 groundModelMatrix;
    XMStoreFloat4x4(&groundModelMatrix, XMMatrixMultiply(
        XMMatrixScaling(1.1f * GRID_WIDTH, 1.f, 1.1f * GRID_DEPTH),
        XMMatrixTranslation(-0.05f * GRID_WIDTH, -1.0f, -0.05f * GRID_DEPTH)
    ));
    modelMatrices.push_back(groundModelMatrix);

    // vector for colors of grid lines
    std::vector<XMFLOAT3> colors = { XMFLOAT3(GROUND_PLANE_COLOR) };

    //cubes for obstacles
    for (const auto& shape : shapes) {
		// if the shape is an obstacle, add it to the scene
        if (shape.functionality == 1) {
            inputStrings.push_back("objs\\cube.obj");
            colors.push_back({ 0.75f, 0.8f, 0.82f });

            XMFLOAT4X4 simShapeMatrix;
            XMStoreFloat4x4(&simShapeMatrix, XMMatrixMultiply(
                XMMatrixScaling(shape.halfSize.x * 2, shape.halfSize.y * 2, shape.halfSize.z * 2),
                XMMatrixTranslation(shape.position.x - shape.halfSize.x, shape.position.y - shape.halfSize.y, shape.position.z - shape.halfSize.z)
            ));
            modelMatrices.push_back(simShapeMatrix);
        }
    }

    //push ground as solid
    auto string = inputStrings.front();
    auto m = modelMatrices.front();
    Mesh newMesh = Mesh((std::filesystem::current_path() / string).string(), context, renderPipeline->getCommandList(), renderPipeline, m, false, colors.front());
    meshes.push_back(newMesh);
    sceneSize += newMesh.getNumTriangles();

    //push shapes as wireframe
    for (int i = 1; i < inputStrings.size(); i++) {
        auto string = inputStrings.at(i);
        auto m = modelMatrices.at(i);
        Mesh newMesh = Mesh((std::filesystem::current_path() / string).string(), context, renderPipeline->getCommandList(), renderPipeline, m, false, colors.at(i));
        meshes.push_back(newMesh);
        sceneSize += newMesh.getNumTriangles();
    }
}

void ObjectScene::draw(Camera* camera) {
    for (Mesh m : meshes) {
        // == IA ==
        auto cmdList = renderPipeline->getCommandList();
        cmdList->IASetVertexBuffers(0, 1, m.getVBV());
        cmdList->IASetIndexBuffer(m.getIBV());

        if (m.getIsWireframe()) {
            cmdList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_LINELIST);
        }
        else {
            cmdList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        }

        // == PSO ==
        cmdList->SetPipelineState(renderPipeline->getPSO());
        cmdList->SetGraphicsRootSignature(renderPipeline->getRootSignature());
        
        // == ROOT ==
        ID3D12DescriptorHeap* descriptorHeaps[] = { renderPipeline->getDescriptorHeap()->GetAddress() };
        cmdList->SetDescriptorHeaps(_countof(descriptorHeaps), descriptorHeaps);
        cmdList->SetGraphicsRootDescriptorTable(1, renderPipeline->getDescriptorHeap()->GetGPUHandleAt(0)); // Descriptor table slot 1 for CBV

        auto viewMat = camera->getViewMat();
        auto projMat = camera->getProjMat();
        cmdList->SetGraphicsRoot32BitConstants(0, 16, &viewMat, 0);
        cmdList->SetGraphicsRoot32BitConstants(0, 16, &projMat, 16);
        cmdList->SetGraphicsRoot32BitConstants(0, 16, m.getModelMatrix(), 32);
        cmdList->SetGraphicsRoot32BitConstants(0, 3, m.getColor(), 48);

        cmdList->DrawIndexedInstanced(m.getNumTriangles() * 3, 1, 0, 0, 0);
    }
}

size_t ObjectScene::getSceneSize() {
    return sceneSize;
}

void ObjectScene::releaseResources() {
	for (Mesh m : meshes) {
		m.releaseResources();
	}
    renderPipeline->releaseResources();
}