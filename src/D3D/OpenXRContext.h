#pragma once
#include "../D3D/DXContext.h"
#include "../Scene/Scene.h"

#define XR_USE_GRAPHICS_API_D3D12
#include <openxr/openxr.h>
#include <openxr/openxr_platform.h>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>

class OpenXRContext {
public:
    struct Swapchain {
        XrSwapchain handle;
        int32_t width;
        int32_t height;
    };

    OpenXRContext();
    ~OpenXRContext();

    XrInstance GetInstance() const { return xrInstance; }
    XrSystemId GetSystemId() const { return systemId; }

    void CreateSession(XrGraphicsBindingD3D12KHR& graphicsBinding);
    XrSession GetSession() const { return xrSession; }

    void CreateSwapchains(ID3D12Device* device);
    void PollEvents(bool& exitRenderLoop, bool& requestRestart);
    void RenderFrame(DXContext* context, ID3D12GraphicsCommandList6* commandList, CommandListID id, Scene& scene, std::vector<CommandListID> commandLists);

private:
    void CreateInstance();

    XrInstance xrInstance = XR_NULL_HANDLE;
    XrSystemId systemId = XR_NULL_SYSTEM_ID;
    XrSession xrSession = XR_NULL_HANDLE;

    uint32_t swapchainWidth = 0;
    uint32_t swapchainHeight = 0;

	std::vector<XrViewConfigurationView> viewConfigViews;
	std::vector<XrView> views;
    int64_t colorSwapchainFormat{ -1 };
    std::vector<Swapchain> swapchains;
    std::map<XrSwapchain, std::vector<XrSwapchainImageBaseHeader*>> swapchainImages;

    std::vector<XrSpace> m_visualizedSpaces;

    ComPointer<ID3D12DescriptorHeap> rtvHeap;
    std::vector<D3D12_CPU_DESCRIPTOR_HANDLE> rtvHandles;

    ComPointer<ID3D12DescriptorHeap> depthStencilHeap;
    ComPointer<ID3D12Resource> depthStencilBuffer;
};
