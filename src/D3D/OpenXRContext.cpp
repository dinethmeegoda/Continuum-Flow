#include "OpenXRContext.h"
#include <iostream>

OpenXRContext::OpenXRContext() {
    CreateInstance();
}

OpenXRContext::~OpenXRContext() {
    if (xrInstance != XR_NULL_HANDLE) {
        xrDestroyInstance(xrInstance);
    }
}

void OpenXRContext::CreateInstance() {
    XrInstanceCreateInfo createInfo{ XR_TYPE_INSTANCE_CREATE_INFO };
    strcpy_s(createInfo.applicationInfo.applicationName, "DX12OpenXREngine");
    createInfo.applicationInfo.apiVersion = XR_CURRENT_API_VERSION;

    std::vector<const char*> extensions = { XR_KHR_D3D12_ENABLE_EXTENSION_NAME };

    createInfo.enabledExtensionCount = static_cast<uint32_t>(extensions.size());
    createInfo.enabledExtensionNames = extensions.data();

    XrResult result = xrCreateInstance(&createInfo, &xrInstance);
    if (XR_FAILED(result)) {
        throw std::runtime_error("Failed to create OpenXR instance.");
    }

    // Select the system (headset)
    XrSystemGetInfo systemInfo{ XR_TYPE_SYSTEM_GET_INFO };
    systemInfo.formFactor = XR_FORM_FACTOR_HEAD_MOUNTED_DISPLAY;

    result = xrGetSystem(xrInstance, &systemInfo, &systemId);
    if (XR_FAILED(result)) {
        throw std::runtime_error("Failed to get OpenXR system.");
    }

    std::cout << "OpenXR Instance Created Successfully!\n";
}

void OpenXRContext::CreateSession(XrGraphicsBindingD3D12KHR& graphicsBinding) {
    // Retrieve the graphics requirements for the OpenXR runtime
    PFN_xrGetD3D12GraphicsRequirementsKHR pfnGetD3D12GraphicsRequirementsKHR = nullptr;
    xrGetInstanceProcAddr(xrInstance, "xrGetD3D12GraphicsRequirementsKHR",
        reinterpret_cast<PFN_xrVoidFunction*>(&pfnGetD3D12GraphicsRequirementsKHR));

    if (!pfnGetD3D12GraphicsRequirementsKHR) {
        throw std::runtime_error("Failed to retrieve xrGetD3D12GraphicsRequirementsKHR function pointer.");
    }

    XrGraphicsRequirementsD3D12KHR graphicsRequirements{ XR_TYPE_GRAPHICS_REQUIREMENTS_D3D12_KHR };
    XrResult result = pfnGetD3D12GraphicsRequirementsKHR(xrInstance, systemId, &graphicsRequirements);
    if (XR_FAILED(result)) {
        throw std::runtime_error("Failed to get OpenXR graphics requirements.");
    }

    // Select a DirectX 12 adapter that matches the OpenXR runtime requirements
    ComPointer<IDXGIFactory1> dxgiFactory;
    if (FAILED(CreateDXGIFactory1(IID_PPV_ARGS(&dxgiFactory)))) {
        throw std::runtime_error("Failed to create DXGI factory.");
    }

    ComPointer<IDXGIAdapter1> dxgiAdapter;
    for (UINT i = 0; dxgiFactory->EnumAdapters1(i, &dxgiAdapter) != DXGI_ERROR_NOT_FOUND; i++) {
        DXGI_ADAPTER_DESC1 adapterDesc;
        dxgiAdapter->GetDesc1(&adapterDesc);

        if (memcmp(&adapterDesc.AdapterLuid, &graphicsRequirements.adapterLuid, sizeof(LUID)) == 0) {
            break; // Found the correct adapter
        }

		dxgiAdapter.Release();
    }

    if (!dxgiAdapter) {
        throw std::runtime_error("Failed to find a compatible DX12 adapter for OpenXR.");
    }

    // Create the OpenXR session
    XrSessionCreateInfo sessionCreateInfo{ XR_TYPE_SESSION_CREATE_INFO };
    sessionCreateInfo.next = &graphicsBinding;
    sessionCreateInfo.systemId = systemId;

    if (XR_FAILED(xrCreateSession(xrInstance, &sessionCreateInfo, &xrSession))) {
        throw std::runtime_error("Failed to create OpenXR session.");
    }

    std::cout << "OpenXR Session Created Successfully!\n";
}

int64_t SelectColorSwapchainFormat(const std::vector<int64_t>& runtimeFormats) {
    // List of supported color swapchain formats.
    constexpr DXGI_FORMAT SupportedColorSwapchainFormats[] = {
        DXGI_FORMAT_R8G8B8A8_UNORM,
        DXGI_FORMAT_B8G8R8A8_UNORM,
        DXGI_FORMAT_R8G8B8A8_UNORM_SRGB,
        DXGI_FORMAT_B8G8R8A8_UNORM_SRGB,
    };

    auto swapchainFormatIt =
        std::find_first_of(runtimeFormats.begin(), runtimeFormats.end(), std::begin(SupportedColorSwapchainFormats),
            std::end(SupportedColorSwapchainFormats));
    if (swapchainFormatIt == runtimeFormats.end()) {
        throw std::runtime_error("No runtime swapchain format supported for color swapchain");
    }

    return *swapchainFormatIt;
}

void OpenXRContext::CreateSwapchains(ID3D12Device* device) {
	// Throw Error if the session is not created
	if (xrSession == XR_NULL_HANDLE) {
		throw std::runtime_error("OpenXR session is not created.");
	}

	// Get the recommended swapchain size
	uint32_t viewCount;
	xrEnumerateViewConfigurationViews(xrInstance, systemId, XR_VIEW_CONFIGURATION_TYPE_PRIMARY_STEREO, 0, &viewCount, nullptr);
	viewConfigViews.resize(viewCount, { XR_TYPE_VIEW_CONFIGURATION_VIEW });
	xrEnumerateViewConfigurationViews(xrInstance, systemId, XR_VIEW_CONFIGURATION_TYPE_PRIMARY_STEREO, viewCount, &viewCount, viewConfigViews.data());

	// Create and cache view buffer for xrLocateViews
	views.resize(viewCount, { XR_TYPE_VIEW });

	// Create the swapchains
    if (viewCount > 0) {
		uint32_t swapchainFormatCount;
		xrEnumerateSwapchainFormats(xrSession, 0, &swapchainFormatCount, nullptr);
		std::vector<int64_t> swapchainFormats(swapchainFormatCount);
        xrEnumerateSwapchainFormats(xrSession, (uint32_t)swapchainFormats.size(), &swapchainFormatCount,
			swapchainFormats.data());
        // Make sure the swapchain format count is the same as the formats size
		if (swapchainFormatCount != swapchainFormats.size()) {
			throw std::runtime_error("Failed to enumerate swapchain formats.");
		}
		colorSwapchainFormat = SelectColorSwapchainFormat(swapchainFormats);

		// Create the swapchains
        for (uint32_t i = 0; i < viewCount; i++) {
			const XrViewConfigurationView& view = viewConfigViews[i];
			XrSwapchainCreateInfo swapchainCreateInfo{ XR_TYPE_SWAPCHAIN_CREATE_INFO };
			swapchainCreateInfo.arraySize = 1;
			swapchainCreateInfo.format = colorSwapchainFormat;
			swapchainCreateInfo.width = view.recommendedImageRectWidth;
			swapchainCreateInfo.height = view.recommendedImageRectHeight;
			swapchainCreateInfo.mipCount = 1;
			swapchainCreateInfo.faceCount = 1;
			swapchainCreateInfo.sampleCount = view.recommendedSwapchainSampleCount;
			swapchainCreateInfo.usageFlags = XR_SWAPCHAIN_USAGE_SAMPLED_BIT | XR_SWAPCHAIN_USAGE_COLOR_ATTACHMENT_BIT;

			Swapchain swapchain;
			swapchain.height = swapchainCreateInfo.height;
			swapchain.width = swapchainCreateInfo.width;
			xrCreateSwapchain(xrSession, &swapchainCreateInfo, &swapchain.handle);

			swapchains.push_back(swapchain);

			uint32_t imageCount;
			xrEnumerateSwapchainImages(swapchain.handle, 0, &imageCount, nullptr);
			std::vector<XrSwapchainImageD3D12KHR> swapchainImages(imageCount, { XR_TYPE_SWAPCHAIN_IMAGE_D3D12_KHR });
        }
    }
}



void OpenXRContext::PollEvents(bool& exitRenderLoop, bool& requestRestart) {
    XrEventDataBuffer eventData{ XR_TYPE_EVENT_DATA_BUFFER };

    while (xrPollEvent(xrInstance, &eventData) == XR_SUCCESS) {
        switch (eventData.type) {
        case XR_TYPE_EVENT_DATA_SESSION_STATE_CHANGED: {
            auto stateEvent = reinterpret_cast<XrEventDataSessionStateChanged*>(&eventData);
            switch (stateEvent->state) {
            case XR_SESSION_STATE_READY: {
                XrSessionBeginInfo beginInfo{ XR_TYPE_SESSION_BEGIN_INFO };
                beginInfo.primaryViewConfigurationType = XR_VIEW_CONFIGURATION_TYPE_PRIMARY_STEREO;
                xrBeginSession(xrSession, &beginInfo);
                break;
            }
            case XR_SESSION_STATE_STOPPING:
                xrEndSession(xrSession);
                break;
            case XR_SESSION_STATE_EXITING:
                exitRenderLoop = true;
                requestRestart = false;
                break;
            case XR_SESSION_STATE_LOSS_PENDING:
                exitRenderLoop = true;
                requestRestart = true;
                break;
            default:
                break;
            }
            break;
        }
        default:
            break;
        }
    }
}

void OpenXRContext::RenderFrame(DXContext* context, ID3D12GraphicsCommandList6* commandList, CommandListID id, Scene& scene, std::vector<CommandListID> commandLists) {
    XrFrameWaitInfo frameWaitInfo{ XR_TYPE_FRAME_WAIT_INFO };
    XrFrameState frameState{ XR_TYPE_FRAME_STATE };
    xrWaitFrame(xrSession, &frameWaitInfo, &frameState);

    XrFrameBeginInfo frameBeginInfo{ XR_TYPE_FRAME_BEGIN_INFO };
    xrBeginFrame(xrSession, &frameBeginInfo);

    std::vector<XrCompositionLayerBaseHeader*> layers;
    XrCompositionLayerProjection layer{ XR_TYPE_COMPOSITION_LAYER_PROJECTION };
    std::vector<XrCompositionLayerProjectionView> projectionViews(swapchains.size());

    for (size_t i = 0; i < swapchains.size(); i++) {
        uint32_t imageIndex;
        XrSwapchainImageAcquireInfo acquireInfo{ XR_TYPE_SWAPCHAIN_IMAGE_ACQUIRE_INFO };
        xrAcquireSwapchainImage(swapchains[i], &acquireInfo, &imageIndex);

        XrSwapchainImageWaitInfo waitInfo{ XR_TYPE_SWAPCHAIN_IMAGE_WAIT_INFO };
        waitInfo.timeout = XR_INFINITE_DURATION;
        xrWaitSwapchainImage(swapchains[i], &waitInfo);

        D3D12_CPU_DESCRIPTOR_HANDLE rtvHandle = rtvHeap->GetCPUDescriptorHandleForHeapStart();
        D3D12_CPU_DESCRIPTOR_HANDLE dsvHandle = depthStencilHeap->GetCPUDescriptorHandleForHeapStart();

        commandList->OMSetRenderTargets(1, &rtvHandle, FALSE, &dsvHandle);

        const float clearColor[4] = { 0.0f, 0.0f, 0.2f, 1.0f };
        commandList->ClearRenderTargetView(rtvHandle, clearColor, 0, nullptr);
        commandList->ClearDepthStencilView(dsvHandle, D3D12_CLEAR_FLAG_DEPTH, 1.0f, 0, 0, nullptr);

        scene.drawPBMPM();
        scene.drawSolidObjects();

        for (auto cmdID : commandLists) {
            context->executeCommandList(cmdID);
        }

        context->executeCommandList(id);

        XrSwapchainImageReleaseInfo releaseInfo{ XR_TYPE_SWAPCHAIN_IMAGE_RELEASE_INFO };
        xrReleaseSwapchainImage(swapchains[i], &releaseInfo);
    }

    layer.viewCount = (uint32_t)projectionViews.size();
    layer.views = projectionViews.data();
    layers.push_back(reinterpret_cast<XrCompositionLayerBaseHeader*>(&layer));

    XrFrameEndInfo frameEndInfo{ XR_TYPE_FRAME_END_INFO };
    frameEndInfo.displayTime = frameState.predictedDisplayTime;
    frameEndInfo.layerCount = (uint32_t)layers.size();
    frameEndInfo.layers = layers.data();

    xrEndFrame(xrSession, &frameEndInfo);
}






