#pragma once
#include "../D3D/DXContext.h"
#include "../Scene/Scene.h"
#include "Scene/SceneConstants.h"

#define XR_USE_GRAPHICS_API_D3D12
#include <openxr/openxr.h>
#include <openxr/openxr_platform.h>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <unordered_map>

#define DEBUG_BREAK __debugbreak()

// XR_DOCS_TAG_BEGIN_Helper_Functions0
inline void OpenXRDebugBreak() {
    std::cerr << "Breakpoint here to debug." << std::endl;
    DEBUG_BREAK;
}

inline const char* GetXRErrorString(XrInstance xrInstance, XrResult result) {
    static char string[XR_MAX_RESULT_STRING_SIZE];
    xrResultToString(xrInstance, result, string);
    return string;
}

inline bool IsStringInVector(std::vector<const char*> list, const char* name) {
    bool found = false;
    for (auto& item : list) {
        if (strcmp(name, item) == 0) {
            found = true;
            break;
        }
    }
    return found;
}

template <typename T>
inline bool BitwiseCheck(const T& value, const T& checkValue) {
    return ((value & checkValue) == checkValue);
}

#define OPENXR_CHECK(x, y)                                                                                                                                  \
    {                                                                                                                                                       \
        XrResult result = (x);                                                                                                                              \
        if (!XR_SUCCEEDED(result)) {                                                                                                                        \
            std::cerr << "ERROR: OPENXR: " << int(result) << "(" << (m_xrInstance ? GetXRErrorString(m_xrInstance, result) : "") << ") " << y << std::endl; \
            OpenXRDebugBreak();                                                                                                                             \
        }                                                                                                                                                   \
    }

#define XR_LOG(...) std::cout << __VA_ARGS__ << "\n"

#define D3D12_CHECK(x, y)                                                                         \
    {                                                                                             \
        HRESULT result = (x);                                                                     \
        if (FAILED(result)) {                                                                     \
            std::cout << "ERROR: D3D12: " << std::hex << "0x" << result << std::dec << std::endl; \
            std::cout << "ERROR: D3D12: " << y << std::endl;                                      \
        }                                                                                         \
    }

#define D3D12_SAFE_RELEASE(p) \
    {                         \
        if (p) {              \
            (p)->Release();   \
            (p) = nullptr;    \
        }                     \
    }

class OpenXRContext {
private:
    enum GraphicsAPI_Type : uint8_t {
        UNKNOWN,
        D3D11,
        D3D12,
        OPENGL,
        OPENGL_ES,
        VULKAN
    };

    struct ImageViewCreateInfo {
        void* image;
        enum class Type : uint8_t {
            RTV,
            DSV,
            SRV,
            UAV
        } type;
        enum class View : uint8_t {
            TYPE_1D,
            TYPE_2D,
            TYPE_3D,
            TYPE_CUBE,
            TYPE_1D_ARRAY,
            TYPE_2D_ARRAY,
            TYPE_CUBE_ARRAY,
        } view;
        int64_t format;
        enum class Aspect : uint8_t {
            COLOR_BIT = 0x01,
            DEPTH_BIT = 0x02,
            STENCIL_BIT = 0x04
        } aspect;
        uint32_t baseMipLevel;
        uint32_t levelCount;
        uint32_t baseArrayLayer;
        uint32_t layerCount;
    };

    struct SwapchainInfo {
        XrSwapchain swapchain = XR_NULL_HANDLE;
        int64_t swapchainFormat = 0;
        std::vector<void*> imageViews;
    };

    enum class SwapchainType : uint8_t {
        COLOR,
        DEPTH
    };

    XrSwapchainImageBaseHeader* AllocateSwapchainImageData(XrSwapchain swapchain, SwapchainType type, uint32_t count);

    void* CreateImageView(const ImageViewCreateInfo& imageViewCI);

    void DestroyImageView(void*& imageView);

    virtual void FreeSwapchainImageData(XrSwapchain swapchain) {
        swapchainImagesMap[swapchain].second.clear();
        swapchainImagesMap.erase(swapchain);
    }

    struct RenderLayerInfo {
        XrTime predictedDisplayTime;
        std::vector<XrCompositionLayerBaseHeader*> layers;
        XrCompositionLayerProjection layerProjection = { XR_TYPE_COMPOSITION_LAYER_PROJECTION };
        std::vector<XrCompositionLayerProjectionView> layerProjectionViews;
    };

    struct Viewport {
        float x;
        float y;
        float width;
        float height;
        float minDepth;
        float maxDepth;
    };
    struct Offset2D {
        int32_t x;
        int32_t y;
    };
    struct Extent2D {
        uint32_t width;
        uint32_t height;
    };
    struct Rect2D {
        Offset2D offset;
        Extent2D extent;
    };

    void SetRenderAttachments(void** colorViews, size_t colorViewCount, void* depthStencilView, uint32_t width, uint32_t height);
    void SetViewports(Viewport* viewports, size_t count);
    void SetScissors(Rect2D* scissors, size_t count);

public:
    struct LeftController {
        float triggerValue = 0.0f;
		float gripValue = 0.0f;
        XMFLOAT3 position = { 0.0f, 0.0f, 0.0f };
        XMFLOAT3 forward = { 0.0f, 0.0f, 1.0f };
    };
    OpenXRContext(ID3D12GraphicsCommandList6* cmdList, DXContext* context, 
        CommandListID id, Camera* c, LeftController &lc);
    ~OpenXRContext();

    void CreateInstance();
    void CreateDebugMessenger();

    void GetInstanceProperties();
    void GetSystemID();

    void CreateSession(XrGraphicsBindingD3D12KHR& graphicsBinding);

    void PollEvents();
    void PollSystemEvents();

	void DestroyDebugMessenger();
	void DestroySession();
    void DestroyInstance();

	bool IsSessionRunning() const { return m_sessionRunning; }
	bool IsApplicationRunning() const { return m_applicationRunning; }

    void GetViewConfigurationViews();
    void CreateSwapchains();
    void DestroySwapchains();

    void GetEnvironmentBlendModes();
    void CreateReferenceSpace();
    void DestroyReferenceSpace();
    void RenderFrame(Scene &scene);
    bool RenderLayer(RenderLayerInfo& renderLayerInfo, Scene &scene);

    void BeginRendering();
	void ClearColor(void* imageView, float r, float g, float b, float a);
	void ClearDepth(void* imageView, float d);
    void EndRendering();

    void CreateActions();

private:
    XrDebugUtilsMessengerEXT CreateOpenXRDebugUtilsMessenger(XrInstance m_xrInstance);
    void DestroyOpenXRDebugUtilsMessenger(XrInstance m_xrInstance, XrDebugUtilsMessengerEXT debugUtilsMessenger);

    XrInstance m_xrInstance = {};
    std::vector<const char*> m_activeAPILayers = {};
    std::vector<const char*> m_activeInstanceExtensions = {};
    std::vector<std::string> m_apiLayers = {};
    std::vector<std::string> m_instanceExtensions = {};

    XrDebugUtilsMessengerEXT m_debugUtilsMessenger = {};

    XrFormFactor m_formFactor = XR_FORM_FACTOR_HEAD_MOUNTED_DISPLAY;
    XrSystemId m_systemID = {};
    XrSystemProperties m_systemProperties = { XR_TYPE_SYSTEM_PROPERTIES };

    XrSession m_session = XR_NULL_HANDLE;

    XrSessionState m_sessionState = XR_SESSION_STATE_UNKNOWN;

    bool m_applicationRunning = true;
    bool m_sessionRunning = false;

    std::vector<XrViewConfigurationType> m_applicationViewConfigurations = { XR_VIEW_CONFIGURATION_TYPE_PRIMARY_STEREO, XR_VIEW_CONFIGURATION_TYPE_PRIMARY_MONO };
    std::vector<XrViewConfigurationType> m_viewConfigurations;
    XrViewConfigurationType m_viewConfiguration = XR_VIEW_CONFIGURATION_TYPE_MAX_ENUM;
    std::vector<XrViewConfigurationView> m_viewConfigurationViews;

    std::vector<SwapchainInfo> m_colorSwapchainInfos = {};
    std::vector<SwapchainInfo> m_depthSwapchainInfos = {};

    std::unordered_map<XrSwapchain, std::pair<SwapchainType, std::vector<XrSwapchainImageD3D12KHR>>> swapchainImagesMap{};

    std::unordered_map<ID3D12Resource*, D3D12_RESOURCE_STATES> imageStates;

    ID3D12Device* m_device = nullptr;

    std::unordered_map<SIZE_T, std::pair<ComPointer<ID3D12DescriptorHeap>, ID3D12Resource*>> imageViewResources;

    virtual void* GetSwapchainImage(XrSwapchain swapchain, uint32_t index) {
        ID3D12Resource* image = swapchainImagesMap[swapchain].second[index].texture;
        D3D12_RESOURCE_STATES state = swapchainImagesMap[swapchain].first == SwapchainType::COLOR ? D3D12_RESOURCE_STATE_RENDER_TARGET : D3D12_RESOURCE_STATE_DEPTH_WRITE;
        imageStates[image] = state;
        return image;
    }

    std::vector<XrEnvironmentBlendMode> m_applicationEnvironmentBlendModes = { XR_ENVIRONMENT_BLEND_MODE_OPAQUE, XR_ENVIRONMENT_BLEND_MODE_ADDITIVE };
    std::vector<XrEnvironmentBlendMode> m_environmentBlendModes = {};
    XrEnvironmentBlendMode m_environmentBlendMode = XR_ENVIRONMENT_BLEND_MODE_MAX_ENUM;

    XrSpace m_localSpace = XR_NULL_HANDLE;

    ID3D12Resource* currentDesktopSwapchainImage = nullptr;

    CommandListID cmdListID;
	DXContext* m_dxContext = nullptr;
    ID3D12GraphicsCommandList6* m_cmdList = nullptr;

	Camera* m_camera = nullptr;

	void UpdateCameraProjectionMatrix(XrView view);

    XrActionSet m_actionSet{};
    XrAction m_moveAction{};
    XrAction m_triggerAction{};
    XrAction m_gripTriggerAction{}; // NEW for middle (grip) trigger
    XrPath m_leftHandPath{};

    XrAction m_leftHandPoseAction{};
    XrSpace m_leftHandSpace{};

    void ApplyCameraMovement(float moveX, float moveZ, float velocity);
    XMFLOAT3 cameraWorldPosition = { 0.0f, 0.7f, 0.0f }; // default player eye height

    LeftController &m_lc;

};
