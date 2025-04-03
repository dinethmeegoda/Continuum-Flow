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

OpenXRContext::OpenXRContext() {
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

    XrSessionCreateInfo sessionCI{ XR_TYPE_SESSION_CREATE_INFO };

    sessionCI.next = &graphicsBinding;
    sessionCI.createFlags = 0;
    sessionCI.systemId = m_systemID;

    OPENXR_CHECK(xrCreateSession(m_xrInstance, &sessionCI, &m_session), "Failed to create Session.");
}

void OpenXRContext::DestroySession() {
    OPENXR_CHECK(xrDestroySession(m_session), "Failed to destroy Session.");
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
                sessionBeginInfo.primaryViewConfigurationType = XR_VIEW_CONFIGURATION_TYPE_PRIMARY_STEREO;
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