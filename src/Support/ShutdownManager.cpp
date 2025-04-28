#include "ShutdownManager.h"

// Static member definitions
std::unique_ptr<DebugLayer> ShutdownManager::debugLayer;
std::unique_ptr<DXContext> ShutdownManager::context;
std::unique_ptr<Scene> ShutdownManager::scene;
std::unique_ptr<Camera> ShutdownManager::camera;
//std::unique_ptr<Keyboard> ShutdownManager::keyboard;
//std::unique_ptr<Mouse> ShutdownManager::mouse;
std::unique_ptr<DescriptorHeap> ShutdownManager::pbmpmDescriptorHeap;
std::unique_ptr<OpenXRContext> ShutdownManager::openxrContext;
std::unique_ptr<OpenXRContext::LeftController> ShutdownManager::leftController;
std::unique_ptr<XrGraphicsBindingD3D12KHR> ShutdownManager::graphicsBinding;

void ShutdownManager::Initialize() {
    debugLayer = std::make_unique<DebugLayer>();
    context = std::make_unique<DXContext>();
    camera = std::make_unique<Camera>();
    //keyboard = std::make_unique<Keyboard>();
    //mouse = std::make_unique<Mouse>();
	pbmpmDescriptorHeap = std::make_unique<DescriptorHeap>(context.get(), D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 150, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE);
    scene = std::make_unique<Scene>(camera.get(), context.get(), OPENXR_CMDLIST_ID, pbmpmDescriptorHeap.get());

    //if (!Window::get().init(context.get(), SCREEN_WIDTH, SCREEN_HEIGHT)) {
    //    //handle could not initialize window
    //    std::cout << "could not initialize window\n";
    //    Window::get().shutdown();
    //}
    leftController = std::make_unique<OpenXRContext::LeftController>();
    graphicsBinding = std::make_unique<XrGraphicsBindingD3D12KHR>();

	openxrContext = std::make_unique<OpenXRContext>(scene.get()->getObjectSolidPipeline()->getCommandList(),
        context.get(), OPENXR_CMDLIST_ID, camera.get(), *leftController.get());
	context->resetCommandList(OPENXR_CMDLIST_ID);

    // Open XR Setup
    openxrContext->CreateInstance();
    openxrContext->CreateDebugMessenger();
    openxrContext->GetInstanceProperties();
    openxrContext->GetSystemID();
    openxrContext->GetViewConfigurationViews();
    openxrContext->GetEnvironmentBlendModes();

    // Create OpenXR session with DX12 graphics binding
    memset(graphicsBinding.get(), 0, sizeof(XrGraphicsBindingD3D12KHR));
    graphicsBinding->type = XR_TYPE_GRAPHICS_BINDING_D3D12_KHR;
    graphicsBinding->device = context->getDevice();
    graphicsBinding->queue = context->getCommandQueue();

	openxrContext->CreateSession(*graphicsBinding.get());
	openxrContext->CreateActions();
    openxrContext->CreateReferenceSpace();
	openxrContext->CreateSwapchains();
}

void ShutdownManager::Shutdown() {
    // Destruct in proper order
    camera.reset();
	//keyboard.reset();
	//mouse.reset();

	// Wait for the GPU to finish
    context->flush(FRAME_COUNT);

	scene.get()->releaseResources();
    scene.reset();

    pbmpmDescriptorHeap->releaseResources();
	pbmpmDescriptorHeap.reset();

    openxrContext->DestroySwapchains();
	openxrContext->DestroyReferenceSpace();
	openxrContext->DestroySession();
	openxrContext->DestroyDebugMessenger();
	openxrContext->DestroyInstance();
	openxrContext.reset();

	graphicsBinding.reset();
	leftController.reset();

    context->flush(1);
    context.reset();
    debugLayer.reset();
}