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