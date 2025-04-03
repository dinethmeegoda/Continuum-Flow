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

class OpenXRContext {
public:

    OpenXRContext();
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
};
