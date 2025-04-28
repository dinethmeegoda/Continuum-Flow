#include "ShutdownManager.h"

// Static member definitions
std::unique_ptr<DebugLayer> ShutdownManager::debugLayer;
std::unique_ptr<DXContext> ShutdownManager::context;
std::unique_ptr<Scene> ShutdownManager::scene;
std::unique_ptr<Camera> ShutdownManager::camera;
std::unique_ptr<Keyboard> ShutdownManager::keyboard;
std::unique_ptr<Mouse> ShutdownManager::mouse;
std::unique_ptr<DescriptorHeap> ShutdownManager::pbmpmDescriptorHeap;

void ShutdownManager::Initialize() {
    debugLayer = std::make_unique<DebugLayer>();
    context = std::make_unique<DXContext>();
    camera = std::make_unique<Camera>();
    keyboard = std::make_unique<Keyboard>();
    mouse = std::make_unique<Mouse>();
	pbmpmDescriptorHeap = std::make_unique<DescriptorHeap>(context.get(), D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 150, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE);
    scene = std::make_unique<Scene>(camera.get(), context.get(), OPENXR_CMDLIST_ID, pbmpmDescriptorHeap.get());

    if (!Window::get().init(context.get(), SCREEN_WIDTH, SCREEN_HEIGHT)) {
        //handle could not initialize window
        std::cout << "could not initialize window\n";
        Window::get().shutdown();
    }
}

void ShutdownManager::Shutdown() {
    // Destruct in proper order
    camera.reset();
	keyboard.reset();
	mouse.reset();

	// Wait for the GPU to finish
    context->flush(FRAME_COUNT);

	scene.get()->releaseResources();
    scene.reset();

    pbmpmDescriptorHeap->releaseResources();
	pbmpmDescriptorHeap.reset();

    context->flush(1);
    context.reset();
    debugLayer.reset();

    Window::get().shutdown();
}