#include "OpenXRContext.h"
#include <iostream>

// Debug Helper Functions

XrBool32 OpenXRMessageCallbackFunction(XrDebugUtilsMessageSeverityFlagsEXT messageSeverity, XrDebugUtilsMessageTypeFlagsEXT messageType, const XrDebugUtilsMessengerCallbackDataEXT* pCallbackData, void* pUserData) {
    // Lambda to covert an XrDebugUtilsMessageSeverityFlagsEXT to std::string. Bitwise check to concatenate multiple severities to the output string.
    auto GetMessageSeverityString = [](XrDebugUtilsMessageSeverityFlagsEXT messageSeverity) -> std::string {
        bool separator = false;

        std::string msgFlags;
        if (BitwiseCheck(messageSeverity, XR_DEBUG_UTILS_MESSAGE_SEVERITY_VERBOSE_BIT_EXT)) {
            msgFlags += "VERBOSE";
            separator = true;
        }
        if (BitwiseCheck(messageSeverity, XR_DEBUG_UTILS_MESSAGE_SEVERITY_INFO_BIT_EXT)) {
            if (separator) {
                msgFlags += ",";
            }
            msgFlags += "INFO";
            separator = true;
        }
        if (BitwiseCheck(messageSeverity, XR_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT)) {
            if (separator) {
                msgFlags += ",";
            }
            msgFlags += "WARN";
            separator = true;
        }
        if (BitwiseCheck(messageSeverity, XR_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT)) {
            if (separator) {
                msgFlags += ",";
            }
            msgFlags += "ERROR";
        }
        return msgFlags;
        };
    // Lambda to covert an XrDebugUtilsMessageTypeFlagsEXT to std::string. Bitwise check to concatenate multiple types to the output string.
    auto GetMessageTypeString = [](XrDebugUtilsMessageTypeFlagsEXT messageType) -> std::string {
        bool separator = false;

        std::string msgFlags;
        if (BitwiseCheck(messageType, XR_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT)) {
            msgFlags += "GEN";
            separator = true;
        }
        if (BitwiseCheck(messageType, XR_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT)) {
            if (separator) {
                msgFlags += ",";
            }
            msgFlags += "SPEC";
            separator = true;
        }
        if (BitwiseCheck(messageType, XR_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT)) {
            if (separator) {
                msgFlags += ",";
            }
            msgFlags += "PERF";
        }
        return msgFlags;
        };

    // Collect message data.
    std::string functionName = (pCallbackData->functionName) ? pCallbackData->functionName : "";
    std::string messageSeverityStr = GetMessageSeverityString(messageSeverity);
    std::string messageTypeStr = GetMessageTypeString(messageType);
    std::string messageId = (pCallbackData->messageId) ? pCallbackData->messageId : "";
    std::string message = (pCallbackData->message) ? pCallbackData->message : "";

    // String stream final message.
    std::stringstream errorMessage;
    errorMessage << functionName << "(" << messageSeverityStr << " / " << messageTypeStr << "): msgNum: " << messageId << " - " << message;

    // Log and debug break.
    std::cerr << errorMessage.str() << std::endl;
    if (BitwiseCheck(messageSeverity, XR_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT)) {
        DEBUG_BREAK;
    }
    return XrBool32();
}

// Swapchain Helper Functions
int64_t SelectColorSwapchainFormat(const std::vector<int64_t>& formats) {
    const std::vector<int64_t>& supportSwapchainFormats = {
        DXGI_FORMAT_R8G8B8A8_UNORM,
        DXGI_FORMAT_B8G8R8A8_UNORM,
        DXGI_FORMAT_R8G8B8A8_UNORM_SRGB,
        DXGI_FORMAT_B8G8R8A8_UNORM_SRGB };

    const std::vector<int64_t>::const_iterator& swapchainFormatIt = std::find_first_of(formats.begin(), formats.end(),
        std::begin(supportSwapchainFormats), std::end(supportSwapchainFormats));
    if (swapchainFormatIt == formats.end()) {
        std::cout << "ERROR: Unable to find supported Color Swapchain Format" << std::endl;
        DEBUG_BREAK;
        return 0;
    }

    return *swapchainFormatIt;
}

int64_t SelectDepthSwapchainFormat(const std::vector<int64_t>& formats) {
    const std::vector<int64_t>& supportSwapchainFormats = {
        DXGI_FORMAT_D32_FLOAT };

    const std::vector<int64_t>::const_iterator& swapchainFormatIt = std::find_first_of(formats.begin(), formats.end(),
        std::begin(supportSwapchainFormats), std::end(supportSwapchainFormats));
    if (swapchainFormatIt == formats.end()) {
        std::cout << "ERROR: Unable to find supported Depth Swapchain Format" << std::endl;
        DEBUG_BREAK;
        return 0;
    }

    return *swapchainFormatIt;
}

XrSwapchainImageBaseHeader* OpenXRContext::AllocateSwapchainImageData(XrSwapchain swapchain, SwapchainType type, uint32_t count) {
    swapchainImagesMap[swapchain].first = type;
    swapchainImagesMap[swapchain].second.resize(count, { XR_TYPE_SWAPCHAIN_IMAGE_D3D12_KHR });
    return reinterpret_cast<XrSwapchainImageBaseHeader*>(swapchainImagesMap[swapchain].second.data());
}

void* OpenXRContext::CreateImageView(const ImageViewCreateInfo& imageViewCI) {
    if (imageViewCI.type == ImageViewCreateInfo::Type::RTV) {
        D3D12_RENDER_TARGET_VIEW_DESC rtvDesc{};
        rtvDesc.Format = (DXGI_FORMAT)imageViewCI.format;

        switch (imageViewCI.view) {
        case ImageViewCreateInfo::View::TYPE_1D: {
            rtvDesc.ViewDimension = D3D12_RTV_DIMENSION_TEXTURE1D;
            rtvDesc.Texture1D.MipSlice = imageViewCI.baseMipLevel;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_2D: {
            rtvDesc.ViewDimension = D3D12_RTV_DIMENSION_TEXTURE2D;
            rtvDesc.Texture2D.MipSlice = imageViewCI.baseMipLevel;
            rtvDesc.Texture2D.PlaneSlice = 0;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_3D: {
            rtvDesc.ViewDimension = D3D12_RTV_DIMENSION_TEXTURE3D;
            rtvDesc.Texture3D.MipSlice = imageViewCI.baseMipLevel;
            rtvDesc.Texture3D.FirstWSlice = imageViewCI.baseArrayLayer;
            rtvDesc.Texture3D.WSize = imageViewCI.layerCount;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_1D_ARRAY: {
            rtvDesc.ViewDimension = D3D12_RTV_DIMENSION_TEXTURE1DARRAY;
            rtvDesc.Texture1DArray.MipSlice = imageViewCI.baseMipLevel;
            rtvDesc.Texture1DArray.FirstArraySlice = imageViewCI.baseArrayLayer;
            rtvDesc.Texture1DArray.ArraySize = imageViewCI.layerCount;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_2D_ARRAY: {
            rtvDesc.ViewDimension = D3D12_RTV_DIMENSION_TEXTURE2DARRAY;
            rtvDesc.Texture2DArray.MipSlice = imageViewCI.baseMipLevel;
            rtvDesc.Texture2DArray.FirstArraySlice = imageViewCI.baseArrayLayer;
            rtvDesc.Texture2DArray.ArraySize = imageViewCI.layerCount;
            rtvDesc.Texture2DArray.PlaneSlice = 0;
            break;
        }
        default:
            DEBUG_BREAK;
            std::cout << "ERROR: D3D12: Unknown ImageView View." << std::endl;
            return nullptr;
        }
        D3D12_CPU_DESCRIPTOR_HANDLE rtv = {};
        ID3D12DescriptorHeap* descHeap;
        D3D12_DESCRIPTOR_HEAP_DESC descHeapDesc;
        descHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV;
        descHeapDesc.NumDescriptors = 1;
        descHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
        descHeapDesc.NodeMask = 0;
        D3D12_CHECK(m_device->CreateDescriptorHeap(&descHeapDesc, IID_PPV_ARGS(&descHeap)), "Failed to create DescriptorHeap.");
        rtv = descHeap->GetCPUDescriptorHandleForHeapStart();
        m_device->CreateRenderTargetView((ID3D12Resource*)imageViewCI.image, &rtvDesc, rtv);
        imageViewResources[rtv.ptr] = { descHeap, (ID3D12Resource*)imageViewCI.image };
        return (void*)rtv.ptr;
    }
    else if (imageViewCI.type == ImageViewCreateInfo::Type::DSV) {
        D3D12_DEPTH_STENCIL_VIEW_DESC dsvDesc{};
        dsvDesc.Format = (DXGI_FORMAT)imageViewCI.format;

        switch (imageViewCI.view) {
        case ImageViewCreateInfo::View::TYPE_1D: {
            dsvDesc.ViewDimension = D3D12_DSV_DIMENSION_TEXTURE1D;
            dsvDesc.Texture1D.MipSlice = imageViewCI.baseMipLevel;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_2D: {
            dsvDesc.ViewDimension = D3D12_DSV_DIMENSION_TEXTURE2D;
            dsvDesc.Texture2D.MipSlice = imageViewCI.baseMipLevel;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_1D_ARRAY: {
            dsvDesc.ViewDimension = D3D12_DSV_DIMENSION_TEXTURE1DARRAY;
            dsvDesc.Texture1DArray.MipSlice = imageViewCI.baseMipLevel;
            dsvDesc.Texture1DArray.FirstArraySlice = imageViewCI.baseArrayLayer;
            dsvDesc.Texture1DArray.ArraySize = imageViewCI.layerCount;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_2D_ARRAY: {
            dsvDesc.ViewDimension = D3D12_DSV_DIMENSION_TEXTURE2DARRAY;
            dsvDesc.Texture2DArray.MipSlice = imageViewCI.baseMipLevel;
            dsvDesc.Texture2DArray.FirstArraySlice = imageViewCI.baseArrayLayer;
            dsvDesc.Texture2DArray.ArraySize = imageViewCI.layerCount;
            break;
        }
        default:
            DEBUG_BREAK;
            std::cout << "ERROR: D3D12: Unknown ImageView View." << std::endl;
            return nullptr;
        }
        D3D12_CPU_DESCRIPTOR_HANDLE dsv = {};
        ID3D12DescriptorHeap* descHeap;
        D3D12_DESCRIPTOR_HEAP_DESC descHeapDesc;
        descHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_DSV;
        descHeapDesc.NumDescriptors = 1;
        descHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
        descHeapDesc.NodeMask = 0;
        D3D12_CHECK(m_device->CreateDescriptorHeap(&descHeapDesc, IID_PPV_ARGS(&descHeap)), "Failed to create DescriptorHeap.");
        dsv = descHeap->GetCPUDescriptorHandleForHeapStart();
        m_device->CreateDepthStencilView((ID3D12Resource*)imageViewCI.image, &dsvDesc, dsv);
        imageViewResources[dsv.ptr] = { descHeap, (ID3D12Resource*)imageViewCI.image };
        return (void*)dsv.ptr;
    }
    else if (imageViewCI.type == ImageViewCreateInfo::Type::SRV) {
        D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc{};
        srvDesc.Format = (DXGI_FORMAT)imageViewCI.format;
        srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;

        switch (imageViewCI.view) {
        case ImageViewCreateInfo::View::TYPE_1D: {
            srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE1D;
            srvDesc.Texture1D.MostDetailedMip = imageViewCI.baseMipLevel;
            srvDesc.Texture1D.MipLevels = imageViewCI.levelCount;
            srvDesc.Texture1D.ResourceMinLODClamp = 0.0f;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_2D: {
            srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2D;
            srvDesc.Texture2D.MostDetailedMip = imageViewCI.baseMipLevel;
            srvDesc.Texture2D.MipLevels = imageViewCI.levelCount;
            srvDesc.Texture2D.PlaneSlice = 0;
            srvDesc.Texture2D.ResourceMinLODClamp = 0.0f;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_3D: {
            srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE3D;
            srvDesc.Texture3D.MostDetailedMip = imageViewCI.baseMipLevel;
            srvDesc.Texture3D.MipLevels = imageViewCI.levelCount;
            srvDesc.Texture3D.ResourceMinLODClamp = 0.0f;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_1D_ARRAY: {
            srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE1DARRAY;
            srvDesc.Texture1DArray.MostDetailedMip = imageViewCI.baseMipLevel;
            srvDesc.Texture1DArray.MipLevels = imageViewCI.levelCount;
            srvDesc.Texture1DArray.FirstArraySlice = imageViewCI.baseArrayLayer;
            srvDesc.Texture1DArray.ArraySize = imageViewCI.layerCount;
            srvDesc.Texture1DArray.ResourceMinLODClamp = 0.0f;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_2D_ARRAY: {
            srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2DARRAY;
            srvDesc.Texture2DArray.MostDetailedMip = imageViewCI.baseMipLevel;
            srvDesc.Texture2DArray.MipLevels = imageViewCI.levelCount;
            srvDesc.Texture2DArray.FirstArraySlice = imageViewCI.baseArrayLayer;
            srvDesc.Texture2DArray.ArraySize = imageViewCI.layerCount;
            srvDesc.Texture2DArray.PlaneSlice = 0;
            srvDesc.Texture2DArray.ResourceMinLODClamp = 0.0f;
            break;
        }
        default:
            DEBUG_BREAK;
            std::cout << "ERROR: D3D12: Unknown ImageView View." << std::endl;
            return nullptr;
        }
        D3D12_CPU_DESCRIPTOR_HANDLE srv = {};
        ID3D12DescriptorHeap* descHeap;
        D3D12_DESCRIPTOR_HEAP_DESC descHeapDesc;
        descHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
        descHeapDesc.NumDescriptors = 1;
        descHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
        descHeapDesc.NodeMask = 0;
        D3D12_CHECK(m_device->CreateDescriptorHeap(&descHeapDesc, IID_PPV_ARGS(&descHeap)), "Failed to create DescriptorHeap.");
        srv = descHeap->GetCPUDescriptorHandleForHeapStart();
        m_device->CreateShaderResourceView((ID3D12Resource*)imageViewCI.image, &srvDesc, srv);
        imageViewResources[srv.ptr] = { descHeap, (ID3D12Resource*)imageViewCI.image };
        return (void*)srv.ptr;
    }
    else if (imageViewCI.type == ImageViewCreateInfo::Type::UAV) {
        D3D12_UNORDERED_ACCESS_VIEW_DESC uavDesc{};
        uavDesc.Format = (DXGI_FORMAT)imageViewCI.format;

        switch (imageViewCI.view) {
        case ImageViewCreateInfo::View::TYPE_1D: {
            uavDesc.ViewDimension = D3D12_UAV_DIMENSION_TEXTURE1D;
            uavDesc.Texture1D.MipSlice = imageViewCI.baseMipLevel;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_2D: {
            uavDesc.ViewDimension = D3D12_UAV_DIMENSION_TEXTURE2D;
            uavDesc.Texture2D.MipSlice = imageViewCI.baseMipLevel;
            uavDesc.Texture2D.PlaneSlice = 0;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_3D: {
            uavDesc.ViewDimension = D3D12_UAV_DIMENSION_TEXTURE3D;
            uavDesc.Texture3D.MipSlice = imageViewCI.baseMipLevel;
            uavDesc.Texture3D.FirstWSlice = imageViewCI.baseArrayLayer;
            uavDesc.Texture3D.WSize = imageViewCI.layerCount;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_1D_ARRAY: {
            uavDesc.ViewDimension = D3D12_UAV_DIMENSION_TEXTURE1DARRAY;
            uavDesc.Texture1DArray.MipSlice = imageViewCI.baseMipLevel;
            uavDesc.Texture1DArray.FirstArraySlice = imageViewCI.baseArrayLayer;
            uavDesc.Texture1DArray.ArraySize = imageViewCI.layerCount;
            break;
        }
        case ImageViewCreateInfo::View::TYPE_2D_ARRAY: {
            uavDesc.ViewDimension = D3D12_UAV_DIMENSION_TEXTURE2DARRAY;
            uavDesc.Texture2DArray.MipSlice = imageViewCI.baseMipLevel;
            uavDesc.Texture2DArray.FirstArraySlice = imageViewCI.baseArrayLayer;
            uavDesc.Texture2DArray.ArraySize = imageViewCI.layerCount;
            uavDesc.Texture2DArray.PlaneSlice = 0;
            break;
        }
        default:
            DEBUG_BREAK;
            std::cout << "ERROR: D3D12: Unknown ImageView View." << std::endl;
            return nullptr;
        }
        D3D12_CPU_DESCRIPTOR_HANDLE uav = {};
        ID3D12DescriptorHeap* descHeap;
        D3D12_DESCRIPTOR_HEAP_DESC descHeapDesc;
        descHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
        descHeapDesc.NumDescriptors = 1;
        descHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
        descHeapDesc.NodeMask = 0;
        D3D12_CHECK(m_device->CreateDescriptorHeap(&descHeapDesc, IID_PPV_ARGS(&descHeap)), "Failed to create DescriptorHeap.");
        uav = descHeap->GetCPUDescriptorHandleForHeapStart();
        m_device->CreateUnorderedAccessView((ID3D12Resource*)imageViewCI.image, nullptr, &uavDesc, uav);
        imageViewResources[uav.ptr] = { descHeap, (ID3D12Resource*)imageViewCI.image };
        return (void*)uav.ptr;
    }
    else {
        DEBUG_BREAK;
        std::cout << "ERROR: D3D12: Unknown ImageView Type." << std::endl;
        return nullptr;
    }
}

void OpenXRContext::DestroyImageView(void*& imageView) {
    D3D12_CPU_DESCRIPTOR_HANDLE d3d12ImageView = { (SIZE_T)imageView };
    ID3D12DescriptorHeap* descHeap = imageViewResources[d3d12ImageView.ptr].first;
    imageViewResources.erase(d3d12ImageView.ptr);
    D3D12_SAFE_RELEASE(descHeap);
    imageView = nullptr;
}

XrDebugUtilsMessengerEXT OpenXRContext::CreateOpenXRDebugUtilsMessenger(XrInstance m_xrInstance) {
    // Fill out a XrDebugUtilsMessengerCreateInfoEXT structure specifying all severities and types.
    // Set the userCallback to OpenXRMessageCallbackFunction().
    XrDebugUtilsMessengerCreateInfoEXT debugUtilsMessengerCI{ XR_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT };
    debugUtilsMessengerCI.messageSeverities = XR_DEBUG_UTILS_MESSAGE_SEVERITY_VERBOSE_BIT_EXT | XR_DEBUG_UTILS_MESSAGE_SEVERITY_INFO_BIT_EXT | XR_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT | XR_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT;
    debugUtilsMessengerCI.messageTypes = XR_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT | XR_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT | XR_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT | XR_DEBUG_UTILS_MESSAGE_TYPE_CONFORMANCE_BIT_EXT;
    debugUtilsMessengerCI.userCallback = (PFN_xrDebugUtilsMessengerCallbackEXT)OpenXRMessageCallbackFunction;
    debugUtilsMessengerCI.userData = nullptr;

    // Load xrCreateDebugUtilsMessengerEXT() function pointer as it is not default loaded by the OpenXR loader.
    XrDebugUtilsMessengerEXT debugUtilsMessenger{};
    PFN_xrCreateDebugUtilsMessengerEXT xrCreateDebugUtilsMessengerEXT;
    OPENXR_CHECK(xrGetInstanceProcAddr(m_xrInstance, "xrCreateDebugUtilsMessengerEXT", (PFN_xrVoidFunction*)&xrCreateDebugUtilsMessengerEXT), "Failed to get InstanceProcAddr.");

    // Finally create and return the XrDebugUtilsMessengerEXT.
    OPENXR_CHECK(xrCreateDebugUtilsMessengerEXT(m_xrInstance, &debugUtilsMessengerCI, &debugUtilsMessenger), "Failed to create DebugUtilsMessenger.");
    return debugUtilsMessenger;
}

void OpenXRContext::DestroyOpenXRDebugUtilsMessenger(XrInstance m_xrInstance, XrDebugUtilsMessengerEXT debugUtilsMessenger) {
    // Load xrDestroyDebugUtilsMessengerEXT() function pointer as it is not default loaded by the OpenXR loader.
    PFN_xrDestroyDebugUtilsMessengerEXT xrDestroyDebugUtilsMessengerEXT;
    OPENXR_CHECK(xrGetInstanceProcAddr(m_xrInstance, "xrDestroyDebugUtilsMessengerEXT", (PFN_xrVoidFunction*)&xrDestroyDebugUtilsMessengerEXT), "Failed to get InstanceProcAddr.");

    // Destroy the provided XrDebugUtilsMessengerEXT.
    OPENXR_CHECK(xrDestroyDebugUtilsMessengerEXT(debugUtilsMessenger), "Failed to destroy DebugUtilsMessenger.");
}

// Rendering Helper Functions
void OpenXRContext::ClearColor(void* imageView, float r, float g, float b, float a) {
    ID3D12Resource* image = imageViewResources[(SIZE_T)imageView].second;
    if (imageStates[image] != D3D12_RESOURCE_STATE_RENDER_TARGET) {
        D3D12_RESOURCE_BARRIER barrier;
        barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
        barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
        barrier.Transition.pResource = image;
        barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
        barrier.Transition.StateBefore = imageStates[currentDesktopSwapchainImage];
        barrier.Transition.StateAfter = imageStates[currentDesktopSwapchainImage] = D3D12_RESOURCE_STATE_RENDER_TARGET;
        m_cmdList->ResourceBarrier(1, &barrier);
    }

    const FLOAT clearColor[4] = { r, g, b, a };
    D3D12_CPU_DESCRIPTOR_HANDLE d3d12ImageView = { (SIZE_T)imageView };
    m_cmdList->ClearRenderTargetView(d3d12ImageView, clearColor, 0, nullptr);
}

void OpenXRContext::ClearDepth(void* imageView, float d) {
    ID3D12Resource* image = imageViewResources[(SIZE_T)imageView].second;
    if (imageStates[image] != D3D12_RESOURCE_STATE_DEPTH_WRITE) {
        D3D12_RESOURCE_BARRIER barrier;
        barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
        barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
        barrier.Transition.pResource = image;
        barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
        barrier.Transition.StateBefore = imageStates[image];
        barrier.Transition.StateAfter = imageStates[image] = D3D12_RESOURCE_STATE_DEPTH_WRITE;
        m_cmdList->ResourceBarrier(1, &barrier);
    }

    D3D12_CPU_DESCRIPTOR_HANDLE d3d12ImageView = { (SIZE_T)imageView };
    m_cmdList->ClearDepthStencilView(d3d12ImageView, D3D12_CLEAR_FLAG_DEPTH, d, 0, 0, nullptr);
}

void OpenXRContext::BeginRendering() {
    //setDescriptorHeap = true;
    //CBV_SRV_UAV_DescriptorOffset = 0;
    //SAMPLER_DescriptorOffset = 0;

    if (currentDesktopSwapchainImage) {
        D3D12_RESOURCE_BARRIER swapchainImageBarrier;
        swapchainImageBarrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
        swapchainImageBarrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
        swapchainImageBarrier.Transition.pResource = currentDesktopSwapchainImage;
        swapchainImageBarrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
        swapchainImageBarrier.Transition.StateBefore = imageStates[currentDesktopSwapchainImage];
        swapchainImageBarrier.Transition.StateAfter = imageStates[currentDesktopSwapchainImage] = D3D12_RESOURCE_STATE_RENDER_TARGET;
        m_cmdList->ResourceBarrier(1, &swapchainImageBarrier);
    }
}

void OpenXRContext::EndRendering() {
    if (currentDesktopSwapchainImage) {
        D3D12_RESOURCE_BARRIER swapchainImageBarrier;
        swapchainImageBarrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
        swapchainImageBarrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
        swapchainImageBarrier.Transition.pResource = currentDesktopSwapchainImage;
        swapchainImageBarrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
        swapchainImageBarrier.Transition.StateBefore = imageStates[currentDesktopSwapchainImage];
        swapchainImageBarrier.Transition.StateAfter = imageStates[currentDesktopSwapchainImage] = D3D12_RESOURCE_STATE_COMMON;
        m_cmdList->ResourceBarrier(1, &swapchainImageBarrier);
    }
    m_dxContext->executeCommandList(cmdListID);
    m_dxContext->resetCommandList(cmdListID);
}

void OpenXRContext::SetRenderAttachments(void** colorViews, size_t colorViewCount,
    void* depthStencilView, uint32_t width, uint32_t height)
{
    std::vector<D3D12_CPU_DESCRIPTOR_HANDLE> d3d12RTVs;
    d3d12RTVs.reserve(colorViewCount);
    for (size_t i = 0; i < colorViewCount; i++) {
        d3d12RTVs.push_back({ (SIZE_T)colorViews[i] });
    }
    D3D12_CPU_DESCRIPTOR_HANDLE d3d12DSV = { (SIZE_T)depthStencilView };

    m_cmdList->OMSetRenderTargets((UINT)colorViewCount, d3d12RTVs.data(), false, &d3d12DSV);
    assert(d3d12DSV.ptr != 0);
}
void OpenXRContext::SetViewports(Viewport* viewports, size_t count) 
{
    std::vector<D3D12_VIEWPORT> d3d12Viewports;
    d3d12Viewports.reserve(count);
    for (size_t i = 0; i < count; i++) {
        const Viewport& viewport = viewports[i];
        d3d12Viewports.push_back({ viewport.x, viewport.y, viewport.width, viewport.height, viewport.minDepth, viewport.maxDepth });
    }

    m_cmdList->RSSetViewports(static_cast<UINT>(d3d12Viewports.size()), d3d12Viewports.data());
}
void OpenXRContext::SetScissors(Rect2D* scissors, size_t count)
{
    std::vector<D3D12_RECT> d3d12Scissors;
    d3d12Scissors.reserve(count);
    for (size_t i = 0; i < count; i++) {
        const Rect2D& scissor = scissors[i];
        d3d12Scissors.push_back({ static_cast<LONG>(scissor.offset.x), static_cast<LONG>(scissor.offset.y), static_cast<LONG>(scissor.extent.width), static_cast<LONG>(scissor.extent.height) });
    }

    m_cmdList->RSSetScissorRects(static_cast<UINT>(d3d12Scissors.size()), d3d12Scissors.data());
}

void OpenXRContext::UpdateCameraProjectionMatrix(XrView headsetView) {
    if (!m_camera) return;

    constexpr float nearZ = 0.05f;
    constexpr float farZ = 100.0f;

    // === PROJECTION MATRIX ===
    XMMATRIX proj;
    {
        const float tanLeft = tanf(headsetView.fov.angleLeft);
        const float tanRight = tanf(headsetView.fov.angleRight);
        const float tanDown = tanf(headsetView.fov.angleDown);
        const float tanUp = tanf(headsetView.fov.angleUp);

        const float left = tanLeft * nearZ;
        const float right = tanRight * nearZ;
        const float bottom = tanDown * nearZ;
        const float top = tanUp * nearZ;

        proj = XMMatrixPerspectiveOffCenterRH(left, right, bottom, top, nearZ, farZ);
    }

    // === VIEW MATRIX ===
    XMMATRIX view;
    {
        const XrVector3f& pos = headsetView.pose.position;
        const XrQuaternionf& rot = headsetView.pose.orientation;

        // Note: OpenXR uses right-handed system, and DirectXMath expects right-handed if we use RH variants
        XMVECTOR headOffset = XMVectorSet(pos.x, pos.y, pos.z, 0.0f);
        XMVECTOR position = XMLoadFloat3(&cameraWorldPosition);
        XMVECTOR worldPos = XMVectorAdd(position, headOffset);
        XMVECTOR orientation = XMVectorSet(rot.x, rot.y, rot.z, rot.w);

        // Transform from local (camera) space to world
        XMMATRIX cameraWorld = XMMatrixRotationQuaternion(orientation) * XMMatrixTranslationFromVector(worldPos);

        // Transform basis vectors from camera space to world space
        XMVECTOR forward = XMVector3Normalize(XMVector3TransformNormal(XMVectorSet(0, 0, -1, 0), cameraWorld));
        XMVECTOR up = XMVector3Normalize(XMVector3TransformNormal(XMVectorSet(0, 1, 0, 0), cameraWorld));
        XMVECTOR right = XMVector3Normalize(XMVector3TransformNormal(XMVectorSet(1, 0, 0, 0), cameraWorld));

        // Store in camera object
        XMStoreFloat3(&m_camera->forward, forward);
        XMStoreFloat3(&m_camera->up, up);
        XMStoreFloat3(&m_camera->right, right);

        // Invert to get view matrix (world -> camera space)
        view = XMMatrixInverse(nullptr, cameraWorld);
    }

    // Store results in camera
    XMStoreFloat4x4(&m_camera->viewMat, view);
    XMStoreFloat4x4(&m_camera->projMat, proj);
}

OpenXRContext::OpenXRContext(ID3D12GraphicsCommandList6* cmdList, 
    DXContext* context, CommandListID commandListID, Camera* camera, LeftController &lc): 
    m_cmdList(cmdList), m_dxContext(context), cmdListID(commandListID), m_camera(camera), m_lc(lc) {
}

OpenXRContext::~OpenXRContext() {
}

void OpenXRContext::CreateInstance() {
    XrApplicationInfo AI;
    strcpy_s(AI.applicationName, "Continuum Flow");
    AI.applicationVersion = 1;
    strcpy_s(AI.engineName, "Breakpoint OpenXR Engine");
    AI.engineVersion = 1;
    AI.apiVersion = XR_CURRENT_API_VERSION;

    m_instanceExtensions.push_back(XR_EXT_DEBUG_UTILS_EXTENSION_NAME);
    m_instanceExtensions.push_back(XR_KHR_D3D12_ENABLE_EXTENSION_NAME);

    // Get all the API Layers from the OpenXR runtime.
    uint32_t apiLayerCount = 0;
    std::vector<XrApiLayerProperties> apiLayerProperties;
    OPENXR_CHECK(xrEnumerateApiLayerProperties(0, &apiLayerCount, nullptr), "Failed to enumerate ApiLayerProperties.");
    apiLayerProperties.resize(apiLayerCount, { XR_TYPE_API_LAYER_PROPERTIES });
    OPENXR_CHECK(xrEnumerateApiLayerProperties(apiLayerCount, &apiLayerCount, apiLayerProperties.data()), "Failed to enumerate ApiLayerProperties.");

    // Check the requested API layers against the ones from the OpenXR. If found add it to the Active API Layers.
    for (auto& requestLayer : m_apiLayers) {
        for (auto& layerProperty : apiLayerProperties) {
            // strcmp returns 0 if the strings match.
            if (strcmp(requestLayer.c_str(), layerProperty.layerName) != 0) {
                continue;
            }
            else {
                m_activeAPILayers.push_back(requestLayer.c_str());
                break;
            }
        }
    }

    // Get all the Instance Extensions from the OpenXR instance.
    uint32_t extensionCount = 0;
    std::vector<XrExtensionProperties> extensionProperties;
    OPENXR_CHECK(xrEnumerateInstanceExtensionProperties(nullptr, 0, &extensionCount, nullptr), "Failed to enumerate InstanceExtensionProperties.");
    extensionProperties.resize(extensionCount, { XR_TYPE_EXTENSION_PROPERTIES });
    OPENXR_CHECK(xrEnumerateInstanceExtensionProperties(nullptr, extensionCount, &extensionCount, extensionProperties.data()), "Failed to enumerate InstanceExtensionProperties.");

    // Check the requested Instance Extensions against the ones from the OpenXR runtime.
    // If an extension is found add it to Active Instance Extensions.
    // Log error if the Instance Extension is not found.
    for (auto& requestedInstanceExtension : m_instanceExtensions) {
        bool found = false;
        for (auto& extensionProperty : extensionProperties) {
            // strcmp returns 0 if the strings match.
            if (strcmp(requestedInstanceExtension.c_str(), extensionProperty.extensionName) != 0) {
                continue;
            }
            else {
                m_activeInstanceExtensions.push_back(requestedInstanceExtension.c_str());
                found = true;
                break;
            }
        }
        if (!found) {
            XR_LOG("Failed to find OpenXR instance extension: " << requestedInstanceExtension);
        }
    }

    XrInstanceCreateInfo instanceCI{ XR_TYPE_INSTANCE_CREATE_INFO };
    instanceCI.createFlags = 0;
    instanceCI.applicationInfo = AI;
    instanceCI.enabledApiLayerCount = static_cast<uint32_t>(m_activeAPILayers.size());
    instanceCI.enabledApiLayerNames = m_activeAPILayers.data();
    instanceCI.enabledExtensionCount = static_cast<uint32_t>(m_activeInstanceExtensions.size());
    instanceCI.enabledExtensionNames = m_activeInstanceExtensions.data();
    OPENXR_CHECK(xrCreateInstance(&instanceCI, &m_xrInstance), "Failed to create Instance.");
}

void OpenXRContext::DestroyInstance() {
    OPENXR_CHECK(xrDestroyInstance(m_xrInstance), "Failed to destroy Instance.");
}

void OpenXRContext::CreateDebugMessenger() {
    // Check that "XR_EXT_debug_utils" is in the active Instance Extensions before creating an XrDebugUtilsMessengerEXT.
    if (IsStringInVector(m_activeInstanceExtensions, XR_EXT_DEBUG_UTILS_EXTENSION_NAME)) {
        m_debugUtilsMessenger = CreateOpenXRDebugUtilsMessenger(m_xrInstance);  // From OpenXRDebugUtils.h.
    }
}

void OpenXRContext::DestroyDebugMessenger() {
    // Check that "XR_EXT_debug_utils" is in the active Instance Extensions before destroying the XrDebugUtilsMessengerEXT.
    if (m_debugUtilsMessenger != XR_NULL_HANDLE) {
        DestroyOpenXRDebugUtilsMessenger(m_xrInstance, m_debugUtilsMessenger);  // From OpenXRDebugUtils.h.
    }
}

void OpenXRContext::GetInstanceProperties() {
    XrInstanceProperties instanceProperties{ XR_TYPE_INSTANCE_PROPERTIES };
    OPENXR_CHECK(xrGetInstanceProperties(m_xrInstance, &instanceProperties), "Failed to get InstanceProperties.");

    XR_LOG("OpenXR Runtime: " << instanceProperties.runtimeName << " - "
        << XR_VERSION_MAJOR(instanceProperties.runtimeVersion) << "."
        << XR_VERSION_MINOR(instanceProperties.runtimeVersion) << "."
        << XR_VERSION_PATCH(instanceProperties.runtimeVersion));
}

void OpenXRContext::GetSystemID() {
    // Get the XrSystemId from the instance and the supplied XrFormFactor.
    XrSystemGetInfo systemGI{ XR_TYPE_SYSTEM_GET_INFO };
    systemGI.formFactor = m_formFactor;
    OPENXR_CHECK(xrGetSystem(m_xrInstance, &systemGI, &m_systemID), "Failed to get SystemID.");

    // Get the System's properties for some general information about the hardware and the vendor.
    OPENXR_CHECK(xrGetSystemProperties(m_xrInstance, m_systemID, &m_systemProperties), "Failed to get SystemProperties.");
}

void OpenXRContext::CreateSession(XrGraphicsBindingD3D12KHR& graphicsBinding) {
    // Retrieve the graphics requirements for the OpenXR runtime
    PFN_xrGetD3D12GraphicsRequirementsKHR pfnGetD3D12GraphicsRequirementsKHR = nullptr;
    xrGetInstanceProcAddr(m_xrInstance, "xrGetD3D12GraphicsRequirementsKHR",
        reinterpret_cast<PFN_xrVoidFunction*>(&pfnGetD3D12GraphicsRequirementsKHR));

    if (!pfnGetD3D12GraphicsRequirementsKHR) {
        throw std::runtime_error("Failed to retrieve xrGetD3D12GraphicsRequirementsKHR function pointer.");
    }

    XrGraphicsRequirementsD3D12KHR graphicsRequirements{ XR_TYPE_GRAPHICS_REQUIREMENTS_D3D12_KHR };
    XrResult result = pfnGetD3D12GraphicsRequirementsKHR(m_xrInstance, m_systemID, &graphicsRequirements);
    if (XR_FAILED(result)) {
        throw std::runtime_error("Failed to get OpenXR graphics requirements.");
    }

	// Record the device in the graphics binding.
	m_device = graphicsBinding.device;

    XrSessionCreateInfo sessionCI{ XR_TYPE_SESSION_CREATE_INFO };

    sessionCI.next = &graphicsBinding;
    sessionCI.createFlags = 0;
    sessionCI.systemId = m_systemID;

    OPENXR_CHECK(xrCreateSession(m_xrInstance, &sessionCI, &m_session), "Failed to create Session.");
}

void OpenXRContext::GetViewConfigurationViews() {
    // Gets the View Configuration Types. The first call gets the count of the array that will be returned. The next call fills out the array.
    uint32_t viewConfigurationCount = 0;
    OPENXR_CHECK(xrEnumerateViewConfigurations(m_xrInstance, m_systemID, 0, &viewConfigurationCount, nullptr), "Failed to enumerate View Configurations.");
    m_viewConfigurations.resize(viewConfigurationCount);
    OPENXR_CHECK(xrEnumerateViewConfigurations(m_xrInstance, m_systemID, viewConfigurationCount, &viewConfigurationCount, m_viewConfigurations.data()), "Failed to enumerate View Configurations.");

    // Pick the first application supported View Configuration Type con supported by the hardware.
    for (const XrViewConfigurationType& viewConfiguration : m_applicationViewConfigurations) {
        if (std::find(m_viewConfigurations.begin(), m_viewConfigurations.end(), viewConfiguration) != m_viewConfigurations.end()) {
            m_viewConfiguration = viewConfiguration;
            break;
        }
    }
    if (m_viewConfiguration == XR_VIEW_CONFIGURATION_TYPE_MAX_ENUM) {
        std::cerr << "Failed to find a view configuration type. Defaulting to XR_VIEW_CONFIGURATION_TYPE_PRIMARY_STEREO." << std::endl;
        m_viewConfiguration = XR_VIEW_CONFIGURATION_TYPE_PRIMARY_STEREO;
    }

    // Gets the View Configuration Views. The first call gets the count of the array that will be returned. The next call fills out the array.
    uint32_t viewConfigurationViewCount = 0;
    OPENXR_CHECK(xrEnumerateViewConfigurationViews(m_xrInstance, m_systemID, m_viewConfiguration, 0, &viewConfigurationViewCount, nullptr), "Failed to enumerate ViewConfiguration Views.");
    m_viewConfigurationViews.resize(viewConfigurationViewCount, { XR_TYPE_VIEW_CONFIGURATION_VIEW });
    OPENXR_CHECK(xrEnumerateViewConfigurationViews(m_xrInstance, m_systemID, m_viewConfiguration, viewConfigurationViewCount, &viewConfigurationViewCount, m_viewConfigurationViews.data()), "Failed to enumerate ViewConfiguration Views.");
}

void OpenXRContext::CreateSwapchains() {
    // Get the supported swapchain formats as an array of int64_t and ordered by runtime preference.
    uint32_t formatCount = 0;
    OPENXR_CHECK(xrEnumerateSwapchainFormats(m_session, 0, &formatCount, nullptr), "Failed to enumerate Swapchain Formats");
    std::vector<int64_t> formats(formatCount);
    OPENXR_CHECK(xrEnumerateSwapchainFormats(m_session, formatCount, &formatCount, formats.data()), "Failed to enumerate Swapchain Formats");
    if (SelectDepthSwapchainFormat(formats) == 0) {
        std::cerr << "Failed to find depth format for Swapchain." << std::endl;
        DEBUG_BREAK;
    }

    //Resize the SwapchainInfo to match the number of view in the View Configuration.
    m_colorSwapchainInfos.resize(m_viewConfigurationViews.size());
    m_depthSwapchainInfos.resize(m_viewConfigurationViews.size());

	// Loop through the View Configuration Views and create a swapchain for each view.
    for (size_t i = 0; i < m_viewConfigurationViews.size(); i++) {
        SwapchainInfo& colorSwapchainInfo = m_colorSwapchainInfos[i];
        SwapchainInfo& depthSwapchainInfo = m_depthSwapchainInfos[i];

        // Fill out an XrSwapchainCreateInfo structure and create an XrSwapchain.
        // Color.
        XrSwapchainCreateInfo swapchainCI{ XR_TYPE_SWAPCHAIN_CREATE_INFO };
        swapchainCI.createFlags = 0;
        swapchainCI.usageFlags = XR_SWAPCHAIN_USAGE_SAMPLED_BIT | XR_SWAPCHAIN_USAGE_COLOR_ATTACHMENT_BIT;
        swapchainCI.format = SelectColorSwapchainFormat(formats);                // Use GraphicsAPI to select the first compatible format.
        swapchainCI.sampleCount = m_viewConfigurationViews[i].recommendedSwapchainSampleCount;  // Use the recommended values from the XrViewConfigurationView.
        swapchainCI.width = m_viewConfigurationViews[i].recommendedImageRectWidth;
        swapchainCI.height = m_viewConfigurationViews[i].recommendedImageRectHeight;
        swapchainCI.faceCount = 1;
        swapchainCI.arraySize = 1;
        swapchainCI.mipCount = 1;
        OPENXR_CHECK(xrCreateSwapchain(m_session, &swapchainCI, &colorSwapchainInfo.swapchain), "Failed to create Color Swapchain");
        colorSwapchainInfo.swapchainFormat = swapchainCI.format;  // Save the swapchain format for later use.

        // Depth.
        swapchainCI.createFlags = 0;
        swapchainCI.usageFlags = XR_SWAPCHAIN_USAGE_SAMPLED_BIT | XR_SWAPCHAIN_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT;
        swapchainCI.format = SelectDepthSwapchainFormat(formats);                // Use GraphicsAPI to select the first compatible format.
        swapchainCI.sampleCount = m_viewConfigurationViews[i].recommendedSwapchainSampleCount;  // Use the recommended values from the XrViewConfigurationView.
        swapchainCI.width = m_viewConfigurationViews[i].recommendedImageRectWidth;
        swapchainCI.height = m_viewConfigurationViews[i].recommendedImageRectHeight;
        swapchainCI.faceCount = 1;
        swapchainCI.arraySize = 1;
        swapchainCI.mipCount = 1;
        OPENXR_CHECK(xrCreateSwapchain(m_session, &swapchainCI, &depthSwapchainInfo.swapchain), "Failed to create Depth Swapchain");
        depthSwapchainInfo.swapchainFormat = swapchainCI.format;  // Save the swapchain format for later use.

        // Get the number of images in the color/depth swapchain and allocate Swapchain image data via GraphicsAPI to store the returned array.
        uint32_t colorSwapchainImageCount = 0;
        OPENXR_CHECK(xrEnumerateSwapchainImages(colorSwapchainInfo.swapchain, 0, &colorSwapchainImageCount, nullptr), "Failed to enumerate Color Swapchain Images.");
        XrSwapchainImageBaseHeader* colorSwapchainImages = AllocateSwapchainImageData(colorSwapchainInfo.swapchain, SwapchainType::COLOR, colorSwapchainImageCount);
        OPENXR_CHECK(xrEnumerateSwapchainImages(colorSwapchainInfo.swapchain, colorSwapchainImageCount, &colorSwapchainImageCount, colorSwapchainImages), "Failed to enumerate Color Swapchain Images.");

        uint32_t depthSwapchainImageCount = 0;
        OPENXR_CHECK(xrEnumerateSwapchainImages(depthSwapchainInfo.swapchain, 0, &depthSwapchainImageCount, nullptr), "Failed to enumerate Depth Swapchain Images.");
        XrSwapchainImageBaseHeader* depthSwapchainImages = AllocateSwapchainImageData(depthSwapchainInfo.swapchain, SwapchainType::DEPTH, depthSwapchainImageCount);
        OPENXR_CHECK(xrEnumerateSwapchainImages(depthSwapchainInfo.swapchain, depthSwapchainImageCount, &depthSwapchainImageCount, depthSwapchainImages), "Failed to enumerate Depth Swapchain Images.");
    
        // Per image in the swapchains, fill out a GraphicsAPI::ImageViewCreateInfo structure and create a color/depth image view.
        for (uint32_t j = 0; j < colorSwapchainImageCount; j++) {
            ImageViewCreateInfo imageViewCI;
            imageViewCI.image = GetSwapchainImage(colorSwapchainInfo.swapchain, j);
            imageViewCI.type = ImageViewCreateInfo::Type::RTV;
            imageViewCI.view = ImageViewCreateInfo::View::TYPE_2D;
            imageViewCI.format = colorSwapchainInfo.swapchainFormat;
            imageViewCI.aspect = ImageViewCreateInfo::Aspect::COLOR_BIT;
            imageViewCI.baseMipLevel = 0;
            imageViewCI.levelCount = 1;
            imageViewCI.baseArrayLayer = 0;
            imageViewCI.layerCount = 1;
            colorSwapchainInfo.imageViews.push_back(CreateImageView(imageViewCI));
        }
        for (uint32_t j = 0; j < depthSwapchainImageCount; j++) {
            ImageViewCreateInfo imageViewCI;
            imageViewCI.image = GetSwapchainImage(depthSwapchainInfo.swapchain, j);
            imageViewCI.type = ImageViewCreateInfo::Type::DSV;
            imageViewCI.view = ImageViewCreateInfo::View::TYPE_2D;
            imageViewCI.format = depthSwapchainInfo.swapchainFormat;
            imageViewCI.aspect = ImageViewCreateInfo::Aspect::DEPTH_BIT;
            imageViewCI.baseMipLevel = 0;
            imageViewCI.levelCount = 1;
            imageViewCI.baseArrayLayer = 0;
            imageViewCI.layerCount = 1;
            depthSwapchainInfo.imageViews.push_back(CreateImageView(imageViewCI));
        }
    }
}

void OpenXRContext::DestroySwapchains() {
    // Per view in the view configuration:
    for (size_t i = 0; i < m_viewConfigurationViews.size(); i++) {
        SwapchainInfo& colorSwapchainInfo = m_colorSwapchainInfos[i];
        SwapchainInfo& depthSwapchainInfo = m_depthSwapchainInfos[i];

        // Destroy the color and depth image views from GraphicsAPI.
        for (void*& imageView : colorSwapchainInfo.imageViews) {
            DestroyImageView(imageView);
        }
        for (void*& imageView : depthSwapchainInfo.imageViews) {
            DestroyImageView(imageView);
        }

        // Free the Swapchain Image Data.
        FreeSwapchainImageData(colorSwapchainInfo.swapchain);
        FreeSwapchainImageData(depthSwapchainInfo.swapchain);

        // Destroy the swapchains.
        OPENXR_CHECK(xrDestroySwapchain(colorSwapchainInfo.swapchain), "Failed to destroy Color Swapchain");
        OPENXR_CHECK(xrDestroySwapchain(depthSwapchainInfo.swapchain), "Failed to destroy Depth Swapchain");
    }
}

void OpenXRContext::DestroySession() {
    OPENXR_CHECK(xrDestroySession(m_session), "Failed to destroy Session.");
}

void OpenXRContext::PollEvents() {
    // Poll OpenXR for a new event.
    XrEventDataBuffer eventData{ XR_TYPE_EVENT_DATA_BUFFER };
    auto XrPollEvents = [&]() -> bool {
        eventData = { XR_TYPE_EVENT_DATA_BUFFER };
        return xrPollEvent(m_xrInstance, &eventData) == XR_SUCCESS;
        };

    while (XrPollEvents()) {
        switch (eventData.type) {
            // Log the number of lost events from the runtime.
        case XR_TYPE_EVENT_DATA_EVENTS_LOST: {
            XrEventDataEventsLost* eventsLost = reinterpret_cast<XrEventDataEventsLost*>(&eventData);
            break;
        }
                                           // Log that an instance loss is pending and shutdown the application.
        case XR_TYPE_EVENT_DATA_INSTANCE_LOSS_PENDING: {
            XrEventDataInstanceLossPending* instanceLossPending = reinterpret_cast<XrEventDataInstanceLossPending*>(&eventData);
            m_sessionRunning = false;
            m_applicationRunning = false;
            break;
        }
                                                     // Log that the interaction profile has changed.
        case XR_TYPE_EVENT_DATA_INTERACTION_PROFILE_CHANGED: {
            XrEventDataInteractionProfileChanged* interactionProfileChanged = reinterpret_cast<XrEventDataInteractionProfileChanged*>(&eventData);
            if (interactionProfileChanged->session != m_session) {
                break;
            }
            break;
        }
                                                           // Log that there's a reference space change pending.
        case XR_TYPE_EVENT_DATA_REFERENCE_SPACE_CHANGE_PENDING: {
            XrEventDataReferenceSpaceChangePending* referenceSpaceChangePending = reinterpret_cast<XrEventDataReferenceSpaceChangePending*>(&eventData);
            if (referenceSpaceChangePending->session != m_session) {
                break;
            }
            break;
        }
                                                              // Session State changes:
        case XR_TYPE_EVENT_DATA_SESSION_STATE_CHANGED: {
            XrEventDataSessionStateChanged* sessionStateChanged = reinterpret_cast<XrEventDataSessionStateChanged*>(&eventData);
            if (sessionStateChanged->session != m_session) {
                break;
            }

            if (sessionStateChanged->state == XR_SESSION_STATE_READY) {
                // SessionState is ready. Begin the XrSession using the XrViewConfigurationType.
                XrSessionBeginInfo sessionBeginInfo{ XR_TYPE_SESSION_BEGIN_INFO };
                sessionBeginInfo.primaryViewConfigurationType = m_viewConfiguration;
                OPENXR_CHECK(xrBeginSession(m_session, &sessionBeginInfo), "Failed to begin Session.");
                m_sessionRunning = true;
            }
            if (sessionStateChanged->state == XR_SESSION_STATE_STOPPING) {
                // SessionState is stopping. End the XrSession.
                OPENXR_CHECK(xrEndSession(m_session), "Failed to end Session.");
                m_sessionRunning = false;
            }
            if (sessionStateChanged->state == XR_SESSION_STATE_EXITING) {
                // SessionState is exiting. Exit the application.
                m_sessionRunning = false;
                m_applicationRunning = false;
            }
            if (sessionStateChanged->state == XR_SESSION_STATE_LOSS_PENDING) {
                // SessionState is loss pending. Exit the application.
                // It's possible to try a reestablish an XrInstance and XrSession, but we will simply exit here
                m_sessionRunning = false;
                m_applicationRunning = false;
            }
            // Store state for reference across the application.
            m_sessionState = sessionStateChanged->state;
            break;
        }
        default: {
            break;
        }
        }
    }
}

void OpenXRContext::PollSystemEvents() {}

void OpenXRContext::GetEnvironmentBlendModes()
{
    // Retrieves the available blend modes. The first call gets the count of the array that will be returned. The next call fills out the array.
    uint32_t environmentBlendModeCount = 0;
    OPENXR_CHECK(xrEnumerateEnvironmentBlendModes(m_xrInstance, m_systemID, m_viewConfiguration, 0, &environmentBlendModeCount, nullptr), "Failed to enumerate EnvironmentBlend Modes.");
    m_environmentBlendModes.resize(environmentBlendModeCount);
    OPENXR_CHECK(xrEnumerateEnvironmentBlendModes(m_xrInstance, m_systemID, m_viewConfiguration, environmentBlendModeCount, &environmentBlendModeCount, m_environmentBlendModes.data()), "Failed to enumerate EnvironmentBlend Modes.");

    // Pick the first application supported blend mode supported by the hardware.
    for (const XrEnvironmentBlendMode& environmentBlendMode : m_applicationEnvironmentBlendModes) {
        if (std::find(m_environmentBlendModes.begin(), m_environmentBlendModes.end(), environmentBlendMode) != m_environmentBlendModes.end()) {
            m_environmentBlendMode = environmentBlendMode;
            break;
        }
    }
    if (m_environmentBlendMode == XR_ENVIRONMENT_BLEND_MODE_MAX_ENUM) {
        XR_LOG("Failed to find a compatible blend mode. Defaulting to XR_ENVIRONMENT_BLEND_MODE_OPAQUE.");
        m_environmentBlendMode = XR_ENVIRONMENT_BLEND_MODE_OPAQUE;
    }
}
void OpenXRContext::CreateReferenceSpace()
{
    // Fill out an XrReferenceSpaceCreateInfo structure and create a reference XrSpace, specifying a Local space with an identity pose as the origin.
    XrReferenceSpaceCreateInfo referenceSpaceCI{ XR_TYPE_REFERENCE_SPACE_CREATE_INFO };
    referenceSpaceCI.referenceSpaceType = XR_REFERENCE_SPACE_TYPE_LOCAL;
    referenceSpaceCI.poseInReferenceSpace = { {0.0f, 0.0f, 0.0f, 1.0f}, {0.0f, 0.0f, 0.0f} };
    OPENXR_CHECK(xrCreateReferenceSpace(m_session, &referenceSpaceCI, &m_localSpace), "Failed to create ReferenceSpace.");
}
void OpenXRContext::DestroyReferenceSpace()
{
    // Destroy the reference XrSpace.
    OPENXR_CHECK(xrDestroySpace(m_localSpace), "Failed to destroy Space.")
}

void OpenXRContext::ApplyCameraMovement(float moveX, float moveZ, float velocity, XrView* headsetView) {
    if (!m_camera) return;

    // Use headset's current rotation to get forward/right directions
    XMMATRIX view = XMLoadFloat4x4(&m_camera->viewMat);
    XMMATRIX world = XMMatrixInverse(nullptr, view);

    XMVECTOR forward = XMVector3Normalize(world.r[2]);  // -Z in world
    XMVECTOR right = XMVector3Normalize(world.r[0]);    // +X in world

    // Remove Y from forward vector so we don't move vertically
    forward = XMVectorSetY(forward, 0.0f);
    forward = XMVector3Normalize(forward);

    right = XMVectorSetY(right, 0.0f);
    right = XMVector3Normalize(right);

    XMVECTOR movement = (-moveZ * forward + moveX * right) * velocity;

    const XrVector3f& pos = headsetView->pose.position;

    // Note: OpenXR uses right-handed system, and DirectXMath expects right-handed if we use RH variants
    XMVECTOR headOffset = XMVectorSet(pos.x, pos.y, pos.z, 0.0f);

	std::cout << "headset pos: " << pos.x << ", " << pos.y << ", " << pos.z << std::endl;

    // Update stored camera position
    XMVECTOR currentPos = XMLoadFloat3(&cameraWorldPosition);
	//currentPos = XMVectorAdd(currentPos, headOffset);
    currentPos = XMVectorAdd(currentPos, movement);
    XMStoreFloat3(&cameraWorldPosition, currentPos);
	m_camera->position = XMFLOAT3(scaleFactorInv * cameraWorldPosition.x, 
                                scaleFactorInv * (cameraWorldPosition.y + pos.y + playerHeight),
                                scaleFactorInv * cameraWorldPosition.z);
}


void OpenXRContext::RenderFrame(Scene &scene)
{
    // Get the XrFrameState for timing and rendering info.
    XrFrameState frameState{ XR_TYPE_FRAME_STATE };
    XrFrameWaitInfo frameWaitInfo{ XR_TYPE_FRAME_WAIT_INFO };
    OPENXR_CHECK(xrWaitFrame(m_session, &frameWaitInfo, &frameState), "Failed to wait for XR Frame.");

    // Tell the OpenXR compositor that the application is beginning the frame.
    XrFrameBeginInfo frameBeginInfo{ XR_TYPE_FRAME_BEGIN_INFO };
    OPENXR_CHECK(xrBeginFrame(m_session, &frameBeginInfo), "Failed to begin the XR Frame.");

    // Variables for rendering and layer composition.
    bool rendered = false;
    RenderLayerInfo renderLayerInfo;
    renderLayerInfo.predictedDisplayTime = frameState.predictedDisplayTime;

    // Check that the session is active and that we should render.
    bool sessionActive = (m_sessionState == XR_SESSION_STATE_SYNCHRONIZED || m_sessionState == XR_SESSION_STATE_VISIBLE || m_sessionState == XR_SESSION_STATE_FOCUSED);
    if (sessionActive && frameState.shouldRender) {
        // Render the stereo image and associate one of swapchain images with the XrCompositionLayerProjection structure.
        rendered = RenderLayer(renderLayerInfo, scene);
        if (rendered) {
            renderLayerInfo.layers.push_back(reinterpret_cast<XrCompositionLayerBaseHeader*>(&renderLayerInfo.layerProjection));
        }
    }

    // Tell OpenXR that we are finished with this frame; specifying its display time, environment blending and layers.
    XrFrameEndInfo frameEndInfo{ XR_TYPE_FRAME_END_INFO };
    frameEndInfo.displayTime = frameState.predictedDisplayTime;
    frameEndInfo.environmentBlendMode = m_environmentBlendMode;
    frameEndInfo.layerCount = static_cast<uint32_t>(renderLayerInfo.layers.size());
    frameEndInfo.layers = renderLayerInfo.layers.data();
    OPENXR_CHECK(xrEndFrame(m_session, &frameEndInfo), "Failed to end the XR Frame.");
}

bool OpenXRContext::RenderLayer(RenderLayerInfo& renderLayerInfo, Scene& scene)
{
    // Locate the views from the view configuration within the (reference) space at the display time.
    std::vector<XrView> views(m_viewConfigurationViews.size(), { XR_TYPE_VIEW });

    XrViewState viewState{ XR_TYPE_VIEW_STATE };  // Will contain information on whether the position and/or orientation is valid and/or tracked.
    XrViewLocateInfo viewLocateInfo{ XR_TYPE_VIEW_LOCATE_INFO };
    viewLocateInfo.viewConfigurationType = m_viewConfiguration;
    viewLocateInfo.displayTime = renderLayerInfo.predictedDisplayTime;
    viewLocateInfo.space = m_localSpace;
    uint32_t viewCount = 0;
    XrResult result = xrLocateViews(m_session, &viewLocateInfo, &viewState, static_cast<uint32_t>(views.size()), &viewCount, views.data());
    if (result != XR_SUCCESS) {
        XR_LOG("Failed to locate Views.");
        return false;
    }

    // Resize the layer projection views to match the view count. The layer projection views are used in the layer projection.
    renderLayerInfo.layerProjectionViews.resize(viewCount, { XR_TYPE_COMPOSITION_LAYER_PROJECTION_VIEW });

    // Per view in the view configuration:
    for (uint32_t i = 0; i < viewCount; i++) {
        SwapchainInfo& colorSwapchainInfo = m_colorSwapchainInfos[i];
        SwapchainInfo& depthSwapchainInfo = m_depthSwapchainInfos[i];

        // Acquire and wait for an image from the swapchains.
        // Get the image index of an image in the swapchains.
        // The timeout is infinite.
        uint32_t colorImageIndex = 0;
        uint32_t depthImageIndex = 0;
        XrSwapchainImageAcquireInfo acquireInfo{ XR_TYPE_SWAPCHAIN_IMAGE_ACQUIRE_INFO };
        OPENXR_CHECK(xrAcquireSwapchainImage(colorSwapchainInfo.swapchain, &acquireInfo, &colorImageIndex), "Failed to acquire Image from the Color Swapchian");
        OPENXR_CHECK(xrAcquireSwapchainImage(depthSwapchainInfo.swapchain, &acquireInfo, &depthImageIndex), "Failed to acquire Image from the Depth Swapchian");

        XrSwapchainImageWaitInfo waitInfo = { XR_TYPE_SWAPCHAIN_IMAGE_WAIT_INFO };
        waitInfo.timeout = XR_INFINITE_DURATION;
        OPENXR_CHECK(xrWaitSwapchainImage(colorSwapchainInfo.swapchain, &waitInfo), "Failed to wait for Image from the Color Swapchain");
        OPENXR_CHECK(xrWaitSwapchainImage(depthSwapchainInfo.swapchain, &waitInfo), "Failed to wait for Image from the Depth Swapchain");

        // Get the width and height and construct the viewport and scissors.
        const uint32_t& width = m_viewConfigurationViews[i].recommendedImageRectWidth;
        const uint32_t& height = m_viewConfigurationViews[i].recommendedImageRectHeight;
        Viewport viewport = { 0.0f, 0.0f, (float)width, (float)height, 0.0f, 1.0f };
        Rect2D scissor = { {(int32_t)0, (int32_t)0}, {width, height} };

        // Fill out the XrCompositionLayerProjectionView structure specifying the pose and fov from the view.
        // This also associates the swapchain image with this layer projection view.
        renderLayerInfo.layerProjectionViews[i] = { XR_TYPE_COMPOSITION_LAYER_PROJECTION_VIEW };
        renderLayerInfo.layerProjectionViews[i].pose = views[i].pose;
        renderLayerInfo.layerProjectionViews[i].fov = views[i].fov;
        renderLayerInfo.layerProjectionViews[i].subImage.swapchain = colorSwapchainInfo.swapchain;
        renderLayerInfo.layerProjectionViews[i].subImage.imageRect.offset.x = 0;
        renderLayerInfo.layerProjectionViews[i].subImage.imageRect.offset.y = 0;
        renderLayerInfo.layerProjectionViews[i].subImage.imageRect.extent.width = static_cast<int32_t>(width);
        renderLayerInfo.layerProjectionViews[i].subImage.imageRect.extent.height = static_cast<int32_t>(height);
        renderLayerInfo.layerProjectionViews[i].subImage.imageArrayIndex = 0;  // Useful for multiview rendering.

        // Rendering code to clear the color and depth image views.
        BeginRendering();

        if (m_environmentBlendMode == XR_ENVIRONMENT_BLEND_MODE_OPAQUE) {
            // VR mode use a background color.
            ClearColor(colorSwapchainInfo.imageViews[colorImageIndex], 0.17f, 0.17f, 0.17f, 1.00f);
        }
        else {
            // In AR mode make the background color black.
            ClearColor(colorSwapchainInfo.imageViews[colorImageIndex], 0.00f, 0.00f, 0.00f, 1.00f);
        }
        ClearDepth(depthSwapchainInfo.imageViews[depthImageIndex], 1.0f);

        // Rendering Stuff
        SetRenderAttachments(&colorSwapchainInfo.imageViews[colorImageIndex], 1, depthSwapchainInfo.imageViews[depthImageIndex], width, height);
        SetViewports(&viewport, 1);
        SetScissors(&scissor, 1);

        // Compute the view-projection transform.
        // All matrices (including OpenXR's) are column-major, right-handed.
        UpdateCameraProjectionMatrix(views[i]);

        // Move
        XrActiveActionSet activeActionSet{};
        activeActionSet.actionSet = m_actionSet;

        XrActionsSyncInfo syncInfo{ XR_TYPE_ACTIONS_SYNC_INFO };
        syncInfo.countActiveActionSets = 1;
        syncInfo.activeActionSets = &activeActionSet;
        OPENXR_CHECK(xrSyncActions(m_session, &syncInfo), "Failed to sync actions");

        // Get the move vector from joystick
        XrActionStateVector2f moveState{ XR_TYPE_ACTION_STATE_VECTOR2F };
        XrActionStateGetInfo getInfo{ XR_TYPE_ACTION_STATE_GET_INFO };
        getInfo.action = m_moveAction;
        OPENXR_CHECK(xrGetActionStateVector2f(m_session, &getInfo, &moveState), "Failed to get move state");

        if (moveState.isActive) {
            float moveX = moveState.currentState.x;
            float moveZ = moveState.currentState.y;

            float deltaTime = 0.01;
            float speedScale = 5;
            // Use moveX and moveZ to update camera/player movement on X and Z axes
            float velocity = std::sqrt(moveState.currentState.x * moveState.currentState.x +
                moveState.currentState.y * moveState.currentState.y) * deltaTime * speedScale;
            ApplyCameraMovement(moveX, moveZ, velocity, &views[i]); // Adjust velocity as needed
        }

        // Check Left Trigger

        XrActionStateFloat triggerState{ XR_TYPE_ACTION_STATE_FLOAT };
        XrActionStateGetInfo getTriggerInfo{ XR_TYPE_ACTION_STATE_GET_INFO };
        getTriggerInfo.action = m_triggerAction;

        OPENXR_CHECK(xrGetActionStateFloat(m_session, &getTriggerInfo, &triggerState), "Failed to get trigger state");

        XrActionStateFloat gripState{ XR_TYPE_ACTION_STATE_FLOAT };
        XrActionStateGetInfo getGripInfo{ XR_TYPE_ACTION_STATE_GET_INFO };
        getGripInfo.action = m_gripTriggerAction;

        OPENXR_CHECK(xrGetActionStateFloat(m_session, &getGripInfo, &gripState), "Failed to get grip trigger state");

        // Consider it pressed if over threshold
        if (triggerState.isActive || gripState.isActive) {

            if (gripState.isActive) {
				m_lc.gripValue = gripState.currentState;
            }
            if (triggerState.isActive) {
                m_lc.triggerValue = triggerState.currentState;
            }

            XrSpaceLocation leftHandLocation{ XR_TYPE_SPACE_LOCATION };
            XrSpaceLocationFlags requiredFlags = XR_SPACE_LOCATION_POSITION_VALID_BIT | XR_SPACE_LOCATION_ORIENTATION_VALID_BIT;

            XrResult result = xrLocateSpace(m_leftHandSpace, m_localSpace, renderLayerInfo.predictedDisplayTime, &leftHandLocation);
            if (XR_SUCCEEDED(result) && (leftHandLocation.locationFlags & requiredFlags) == requiredFlags) {
                const XrPosef& pose = leftHandLocation.pose;

                // World-space position of the controller
                m_lc.position = { (pose.position.x + cameraWorldPosition.x) * scaleFactorInv, (pose.position.y + cameraWorldPosition.y + playerHeight) * scaleFactorInv, (pose.position.z + cameraWorldPosition.z) * scaleFactorInv };
                
                // Orientation -> direction (as covered earlier)
                XMVECTOR orientation = XMVectorSet(pose.orientation.x, pose.orientation.y, pose.orientation.z, pose.orientation.w);
                XMVECTOR localForward = XMVectorSet(0, 0, -1, 0);
                XMVECTOR worldForward = XMVector3Rotate(localForward, orientation);
                worldForward = XMVector3Normalize(worldForward);

                XMStoreFloat3(&m_lc.forward, worldForward);
            }
        }

        scene.drawSolidObjects();
        //scene.drawSpawners();
        scene.drawPBMPM();
        scene.drawFluid(0, 0);

        EndRendering();

        // Give the swapchain image back to OpenXR, allowing the compositor to use the image.
        XrSwapchainImageReleaseInfo releaseInfo{ XR_TYPE_SWAPCHAIN_IMAGE_RELEASE_INFO };
        OPENXR_CHECK(xrReleaseSwapchainImage(colorSwapchainInfo.swapchain, &releaseInfo), "Failed to release Image back to the Color Swapchain");
        OPENXR_CHECK(xrReleaseSwapchainImage(depthSwapchainInfo.swapchain, &releaseInfo), "Failed to release Image back to the Depth Swapchain");
    }

    // Fill out the XrCompositionLayerProjection structure for usage with xrEndFrame().
    renderLayerInfo.layerProjection.layerFlags = XR_COMPOSITION_LAYER_BLEND_TEXTURE_SOURCE_ALPHA_BIT | XR_COMPOSITION_LAYER_CORRECT_CHROMATIC_ABERRATION_BIT;
    renderLayerInfo.layerProjection.space = m_localSpace;
    renderLayerInfo.layerProjection.viewCount = static_cast<uint32_t>(renderLayerInfo.layerProjectionViews.size());
    renderLayerInfo.layerProjection.views = renderLayerInfo.layerProjectionViews.data();

    return true;
}

void OpenXRContext::CreateActions() {
    // === Create Action Set ===
    XrActionSetCreateInfo actionSetInfo{ XR_TYPE_ACTION_SET_CREATE_INFO };
    strcpy_s(actionSetInfo.actionSetName, "main_action_set");
    strcpy_s(actionSetInfo.localizedActionSetName, "Main Action Set");
    actionSetInfo.priority = 0;
    OPENXR_CHECK(xrCreateActionSet(m_xrInstance, &actionSetInfo, &m_actionSet), "Failed to create action set");

    // === Define Subaction Path for Left Hand ===
    OPENXR_CHECK(xrStringToPath(m_xrInstance, "/user/hand/left", &m_leftHandPath), "Failed to get left hand path");

    // === Movement Action (Vector2f) ===
    XrActionCreateInfo moveActionInfo{ XR_TYPE_ACTION_CREATE_INFO };
    moveActionInfo.actionType = XR_ACTION_TYPE_VECTOR2F_INPUT;
    strcpy_s(moveActionInfo.actionName, "move");
    strcpy_s(moveActionInfo.localizedActionName, "Move");
    moveActionInfo.countSubactionPaths = 1;
    moveActionInfo.subactionPaths = &m_leftHandPath;
    OPENXR_CHECK(xrCreateAction(m_actionSet, &moveActionInfo, &m_moveAction), "Failed to create move action");

    // === Index Trigger Action (Float Input) ===
    XrActionCreateInfo triggerActionInfo{ XR_TYPE_ACTION_CREATE_INFO };
    triggerActionInfo.actionType = XR_ACTION_TYPE_FLOAT_INPUT;
    strcpy_s(triggerActionInfo.actionName, "left_trigger");
    strcpy_s(triggerActionInfo.localizedActionName, "Left Trigger");
    triggerActionInfo.countSubactionPaths = 1;
    triggerActionInfo.subactionPaths = &m_leftHandPath;
    OPENXR_CHECK(xrCreateAction(m_actionSet, &triggerActionInfo, &m_triggerAction), "Failed to create trigger action");

    // === Middle Trigger (Grip) Action (Float Input) ===
    XrActionCreateInfo gripActionInfo{ XR_TYPE_ACTION_CREATE_INFO };
    gripActionInfo.actionType = XR_ACTION_TYPE_FLOAT_INPUT;
    strcpy_s(gripActionInfo.actionName, "left_grip");
    strcpy_s(gripActionInfo.localizedActionName, "Left Grip");
    gripActionInfo.countSubactionPaths = 1;
    gripActionInfo.subactionPaths = &m_leftHandPath;
    OPENXR_CHECK(xrCreateAction(m_actionSet, &gripActionInfo, &m_gripTriggerAction), "Failed to create grip action");

    // === Pose Action (for controller tracking) ===
    XrActionCreateInfo poseActionInfo{ XR_TYPE_ACTION_CREATE_INFO };
    poseActionInfo.actionType = XR_ACTION_TYPE_POSE_INPUT;
    strcpy_s(poseActionInfo.actionName, "left_hand_pose");
    strcpy_s(poseActionInfo.localizedActionName, "Left Hand Pose");
    poseActionInfo.countSubactionPaths = 1;
    poseActionInfo.subactionPaths = &m_leftHandPath;
    OPENXR_CHECK(xrCreateAction(m_actionSet, &poseActionInfo, &m_leftHandPoseAction), "Failed to create left hand pose action");

    // === Suggest Bindings ===
    XrPath thumbstickPath, triggerValuePath, gripValuePath, gripPosePath;
    OPENXR_CHECK(xrStringToPath(m_xrInstance, "/user/hand/left/input/thumbstick", &thumbstickPath), "Failed to get thumbstick path");
    OPENXR_CHECK(xrStringToPath(m_xrInstance, "/user/hand/left/input/trigger/value", &triggerValuePath), "Failed to get trigger path");
    OPENXR_CHECK(xrStringToPath(m_xrInstance, "/user/hand/left/input/squeeze/value", &gripValuePath), "Failed to get squeeze (grip) path");
    OPENXR_CHECK(xrStringToPath(m_xrInstance, "/user/hand/left/input/grip/pose", &gripPosePath), "Failed to get grip pose path");

    XrActionSuggestedBinding bindings[] = {
        { m_moveAction, thumbstickPath },
        { m_triggerAction, triggerValuePath },
        { m_gripTriggerAction, gripValuePath },
        { m_leftHandPoseAction, gripPosePath }
    };

    XrInteractionProfileSuggestedBinding suggestedBindings{ XR_TYPE_INTERACTION_PROFILE_SUGGESTED_BINDING };
    OPENXR_CHECK(xrStringToPath(m_xrInstance, "/interaction_profiles/oculus/touch_controller", &suggestedBindings.interactionProfile), "Failed to get interaction profile path");
    suggestedBindings.suggestedBindings = bindings;
    suggestedBindings.countSuggestedBindings = (uint32_t)std::size(bindings);

    OPENXR_CHECK(xrSuggestInteractionProfileBindings(m_xrInstance, &suggestedBindings), "Failed to suggest bindings");

    // === Attach Action Set to Session ===
    XrSessionActionSetsAttachInfo attachInfo{ XR_TYPE_SESSION_ACTION_SETS_ATTACH_INFO };
    attachInfo.countActionSets = 1;
    attachInfo.actionSets = &m_actionSet;
    OPENXR_CHECK(xrAttachSessionActionSets(m_session, &attachInfo), "Failed to attach action set");

    // === Create Action Space for Pose Tracking ===
    XrActionSpaceCreateInfo spaceCreateInfo{ XR_TYPE_ACTION_SPACE_CREATE_INFO };
    spaceCreateInfo.action = m_leftHandPoseAction;
    spaceCreateInfo.subactionPath = m_leftHandPath;
    spaceCreateInfo.poseInActionSpace = { {0,0,0,1}, {0,0,0} }; // Identity pose
    OPENXR_CHECK(xrCreateActionSpace(m_session, &spaceCreateInfo, &m_leftHandSpace), "Failed to create left hand space");
}



