#pragma once

#include "Scene/Camera.h"

#include "D3D/Pipeline/RenderPipeline.h"
#include "D3D/Pipeline/MeshPipeline.h"

class Drawable {
public:
	Drawable() = delete;
	Drawable(DXContext* context, RenderPipeline* pipeline);
	Drawable(DXContext* context, MeshPipeline* pipeline);

	void constructScene();

	void draw(Camera* camera);

	void releaseResources();
	const float playerHeight = 1.6f;
	const float scaleFactor = 0.1f;

protected:
	DXContext* context;
	RenderPipeline* renderPipeline;
	MeshPipeline* meshPipeline;
};