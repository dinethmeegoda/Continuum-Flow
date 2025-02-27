#pragma once
#include "../Support/WinInclude.h"
#include "../Support/ComPointer.h"

#define XR_USE_GRAPHICS_API_D3D12
#include <openxr/openxr.h>
#include <openxr/openxr_platform.h>
#include <vector>
#include <string>
#include <stdexcept>

class OpenXRContext {
public:
    OpenXRContext();
    ~OpenXRContext();

    XrInstance GetInstance() const { return xrInstance; }
    XrSystemId GetSystemId() const { return systemId; }

    void CreateSession(XrGraphicsBindingD3D12KHR& graphicsBinding);
    XrSession GetSession() const { return xrSession; }

private:
    void CreateInstance();

    XrInstance xrInstance = XR_NULL_HANDLE;
    XrSystemId systemId = XR_NULL_SYSTEM_ID;
    XrSession xrSession = XR_NULL_HANDLE;
};
