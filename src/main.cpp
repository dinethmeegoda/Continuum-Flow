#include "main.h"

int main() {
    //initialize scene
	ShutdownManager::Initialize();
	// Get Pointers to the context, camera, global descriptor heap, scene, and XR components
	DXContext* context = ShutdownManager::GetContext();
	Camera* camera = ShutdownManager::GetCamera();
	DescriptorHeap* pbmpmDescriptorHeap = ShutdownManager::GetPBMPMDescriptorHeap();
	Scene* scene = ShutdownManager::GetScene();
	OpenXRContext* openXR = ShutdownManager::GetOpenXRContext();
	OpenXRContext::Controller* leftController = ShutdownManager::GetLeftController();
	OpenXRContext::Controller* rightController = ShutdownManager::GetRightController();

    PBMPMConstants* pbmpmCurrConstants = scene->getPBMPMConstants();
	InteractionConstants* pbmpmInteractionConstants = scene->getPBMPMInteractionConstants();
    float maxForceStrength = 10.0f;

    unsigned int renderOptions = 0;

    while (openXR->IsApplicationRunning()) {

		    // Check for enable gravity
		     if (openXR->enableGravity) {
			     pbmpmCurrConstants->gravityStrength = 2.5f;
		     }
             else {
				 pbmpmCurrConstants->gravityStrength = 0.0f;
			 }
             if (leftController->triggerValue > 0.05 || leftController->gripValue > 0.05) {
                 // enable interaction force
                 pbmpmInteractionConstants->leftActivation = 1;

                 // Pulling Fluid
                 if (leftController->triggerValue > 0.05) {
                     pbmpmInteractionConstants->leftFunction = 0;
                     pbmpmInteractionConstants->leftStrength = maxForceStrength * leftController->triggerValue;
                     //std::cout << "Trigger Value: " << pbmpmIterConstants.mouseStrength << std::endl;
                 }
                 if (leftController->gripValue > 0.05) {
                     pbmpmInteractionConstants->leftFunction = 1;
                     pbmpmInteractionConstants->leftStrength = maxForceStrength * leftController->gripValue;
                     //std::cout << "Trigger Value: " << pbmpmIterConstants.mouseStrength << std::endl;
                 }

                 pbmpmInteractionConstants->leftPosition = XMFLOAT3(leftController->position.x * scaleFactorInv,
                     (leftController->position.y + playerHeight) * scaleFactorInv, leftController->position.z * scaleFactorInv);
                 pbmpmInteractionConstants->leftRayDirection = XMFLOAT3(leftController->forward.x,
                     leftController->forward.y, leftController->forward.z);
             }
             else {
                 pbmpmInteractionConstants->leftActivation = 0;
             }

             if (rightController->triggerValue > 0.05 || rightController->gripValue > 0.05) {
                 // enable interaction force
                 pbmpmInteractionConstants->rightActivation = 1;

                 // Pulling Fluid
                 if (rightController->triggerValue > 0.05) {
                     pbmpmInteractionConstants->rightFunction = 0;
                     pbmpmInteractionConstants->rightStrength = maxForceStrength * rightController->triggerValue;
                 }
                 if (rightController->gripValue > 0.05) {
                     pbmpmInteractionConstants->rightFunction = 1;
                     pbmpmInteractionConstants->rightStrength = maxForceStrength * rightController->gripValue;
                 }

                 pbmpmInteractionConstants->rightPosition = XMFLOAT3(rightController->position.x * scaleFactorInv,
                     (rightController->position.y + playerHeight) * scaleFactorInv, rightController->position.z * scaleFactorInv);
                 pbmpmInteractionConstants->rightRayDirection = XMFLOAT3(rightController->forward.x,
                     rightController->forward.y, rightController->forward.z);
             }
             else {
                 pbmpmInteractionConstants->rightActivation = 0;
             }

             //compute pbmpm + mesh shader
             int renderMode = openXR->getRenderMode();
             scene->compute(renderMode == 0 || renderMode == 1);

        	 openXR->PollSystemEvents();
             openXR->PollEvents();

             if (openXR->IsSessionRunning()) {
                 // Render Frame
                 openXR->RenderFrame(*scene);
             }
    }

    // Now you can safely call ReportLiveDeviceObjects
   /* ComPointer<ID3D12DebugDevice> debugDevice;
    if (SUCCEEDED(context->getDevice()->QueryInterface(IID_PPV_ARGS(&debugDevice)))) {
    }*/

	ShutdownManager::Shutdown();
    return 0;
    /*if (debugDevice) {
        debugDevice->ReportLiveDeviceObjects(D3D12_RLDO_DETAIL);
    }*/

    /*ImGui_ImplDX12_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImGui::DestroyContext();
    imguiSRVHeap->Release();*/

}
