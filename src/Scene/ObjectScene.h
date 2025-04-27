#pragma once

#include <filesystem>
#include <vector>

#include "Drawable.h"
#include "Mesh.h"
#include "PBMPMScene.h"

class ObjectScene : public Drawable {
public:
	ObjectScene(DXContext* context, RenderPipeline* pipeline, std::vector<SimShape>& shapes, int renderWireframe = 0);

	void constructSceneGrid();
	void constructSceneSpawners();
	void constructSceneSolid();
	void constructLasers();

	void draw(Camera* camera);
	void drawLasers(Camera* camera, XMMATRIX* leftLaser);

	size_t getSceneSize();

	void releaseResources();

	bool instanced{ false };

private:
	std::vector<Mesh> meshes;
	std::vector<XMFLOAT4X4> modelMatrices;

	std::vector<SimShape> shapes;

	size_t sceneSize{ 0 };
};