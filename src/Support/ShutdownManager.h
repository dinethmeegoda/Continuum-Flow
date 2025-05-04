#pragma once

#include <memory>
#include "Debug/DebugLayer.h"
#include "D3D/DXContext.h"
#include "D3D/DescriptorHeap.h"
#include "Scene/Scene.h"
#include "Window.h"
#include "D3D/OpenXRContext.h"

class ShutdownManager {
public:
    static void Initialize();
    static void Shutdown();

	static DebugLayer* GetDebugLayer() { return debugLayer.get(); }
	static DXContext* GetContext() { return context.get(); }
	static Scene* GetScene() { return scene.get(); }
	static Camera* GetCamera() { return camera.get(); }
	//static Keyboard* GetKeyboard() { return keyboard.get(); }
	//static Mouse* GetMouse() { return mouse.get(); }
	static DescriptorHeap* GetPBMPMDescriptorHeap() { return pbmpmDescriptorHeap.get(); }
	static OpenXRContext* GetOpenXRContext() { return openxrContext.get(); }
	static OpenXRContext::Controller* GetLeftController() { return leftController.get(); }
	static OpenXRContext::Controller* GetRightController() { return rightController.get(); }

private:
    static std::unique_ptr<DebugLayer> debugLayer;
    static std::unique_ptr<DXContext> context;
    static std::unique_ptr<Scene> scene;
    static std::unique_ptr<Camera> camera;
	static std::unique_ptr<OpenXRContext> openxrContext;
    //static std::unique_ptr<Keyboard> keyboard;
    //static std::unique_ptr<Mouse> mouse;
	static std::unique_ptr<DescriptorHeap> pbmpmDescriptorHeap;
	static std::unique_ptr<OpenXRContext::Controller> leftController;
	static std::unique_ptr<OpenXRContext::Controller> rightController;
	static std::unique_ptr<XrGraphicsBindingD3D12KHR> graphicsBinding;
};
