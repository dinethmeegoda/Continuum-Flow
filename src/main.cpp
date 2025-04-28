#include "main.h"

int main() {
    //set up DX, window, keyboard mouse
    // Create global descriptor heap
    //initialize scene
	ShutdownManager::Initialize();
	DXContext* context = ShutdownManager::GetContext();
	Camera* camera = ShutdownManager::GetCamera();
	//Keyboard* keyboard = ShutdownManager::GetKeyboard();
	//Mouse* mouse = ShutdownManager::GetMouse();
	DescriptorHeap* pbmpmDescriptorHeap = ShutdownManager::GetPBMPMDescriptorHeap();
	Scene* scene = ShutdownManager::GetScene();
	OpenXRContext* openXR = ShutdownManager::GetOpenXRContext();
	OpenXRContext::LeftController* leftController = ShutdownManager::GetLeftController();

    // Create Left Controller Data Struct
    //OpenXRContext::LeftController leftController;

    // Initialize OpenXR
    /*OpenXRContext openXR(scene.getObjectSolidPipeline()->getCommandList(),
        &context, OPENXR_CMDLIST_ID, camera.get(), leftController);
    context.resetCommandList(OPENXR_CMDLIST_ID);*/

    // Create OpenXR instance
    //openXR.CreateInstance();

    // Create OpenXR Debug Messager and log Instance Properties & System ID
    //openXR.CreateDebugMessenger();
    //openXR.GetInstanceProperties();
    //openXR.GetSystemID();

    //openXR.GetViewConfigurationViews();
    //openXR.GetEnvironmentBlendModes();

    // Create OpenXR session with DX12 graphics binding
    //XrGraphicsBindingD3D12KHR graphicsBinding{ XR_TYPE_GRAPHICS_BINDING_D3D12_KHR };
    //graphicsBinding.device = context.getDevice();
    //graphicsBinding.queue = context.getCommandQueue();

    //openXR.CreateSession(graphicsBinding);
    //openXR.CreateActions();

    //std::cout << "DX12 Engine with OpenXR Initialized Successfully!\n";

    // Create OpenXR reference space
    //openXR.CreateReferenceSpace();

    // Create Swapchains for OpenXR
    //openXR.CreateSwapchains();

    //initialize ImGUI
    //ImGuiIO& io = initImGUI(*context);

    //set mouse to use the window
    //mouse->SetWindow(Window::get().getHWND());

    // Get the client area of the window
    /*RECT rect;
    GetClientRect(Window::get().getHWND(), &rect);
    float clientWidth = static_cast<float>(rect.right - rect.left);
    float clientHeight = static_cast<float>(rect.bottom - rect.top);*/

    PBMPMConstants pbmpmCurrConstants = scene->getPBMPMConstants();
    PBMPMConstants pbmpmIterConstants = pbmpmCurrConstants;
    float maxForceStrength = 10.0f;

    unsigned int renderOptions = 0;

    //bool exitRenderLoop = false, requestRestart = false;
    //

    while (openXR->IsApplicationRunning()) {
        //update window
        /*Window::get().update();
        if (Window::get().getShouldResize()) {
            //flush pending buffer operations in swapchain
            context->flush(FRAME_COUNT);
            Window::get().resize();
            camera->updateAspect((float)Window::get().getWidth() / (float)Window::get().getHeight());
        }

        auto kState = keyboard->GetState();
        auto mState = mouse->GetState();
        mouse->SetMode(mState.leftButton ? Mouse::MODE_RELATIVE : Mouse::MODE_ABSOLUTE);
        camera->kmStateCheck(kState, mState);

        if (mState.rightButton) {

            // If right mouse button is pressed, we should update constants

            if (kState.LeftShift) {
                // Pulling Fluid
                pbmpmIterConstants.mouseFunction = 2;
            }
            else if (kState.LeftAlt) {
                // Grab Fluid Ball
                pbmpmIterConstants.mouseFunction = 1;
            }
            else {
                // Pushing Fluid
                pbmpmIterConstants.mouseFunction = 0;
            }

            // enable mouse force
            pbmpmIterConstants.mouseActivation = 1;

            POINT mousePos;
            GetCursorPos(&mousePos);
            ScreenToClient(Window::get().getHWND(), &mousePos);
            float ndcX = (2.0f * mousePos.x / clientWidth) - 1.0f;
            float ndcY = 1.0f - (2.0f * mousePos.y / clientHeight);

            XMFLOAT4 prevMousePos = pbmpmIterConstants.mousePosition;
            ComputeMouseRay(
                Window::get().getHWND(),
                ndcX,
                ndcY,
                camera->getProjMat(),
                camera->getViewMat(),
                pbmpmIterConstants.mousePosition,
                pbmpmIterConstants.mouseRayDirection
            );
        }
        else {
            pbmpmIterConstants.mouseActivation = 0;
        }*/

             if (leftController->triggerValue > 0.05 || leftController->gripValue > 0.05) {
                 // enable mouse force
                 pbmpmIterConstants.mouseActivation = 1;

                 // Pulling Fluid
                 if (leftController->triggerValue > 0.05) {
                     pbmpmIterConstants.mouseFunction = 0;
                 	pbmpmIterConstants.mouseStrength = maxForceStrength * leftController->triggerValue;
                 	//std::cout << "Trigger Value: " << pbmpmIterConstants.mouseStrength << std::endl;
                 }
                 if (leftController->gripValue > 0.05) {
                 	pbmpmIterConstants.mouseFunction = 1;
                     pbmpmIterConstants.mouseStrength = maxForceStrength * leftController->gripValue;
                     //std::cout << "Trigger Value: " << pbmpmIterConstants.mouseStrength << std::endl;
                 }

                 pbmpmIterConstants.mousePosition = XMFLOAT4(leftController->position.x,
                     leftController->position.y, leftController->position.z, 1.0);
                 pbmpmIterConstants.mouseRayDirection = XMFLOAT4(leftController->forward.x,
                     leftController->forward.y, leftController->forward.z, 1.0);
             }
             else {
                 pbmpmIterConstants.mouseActivation = 0;
             }

             //compute pbmpm + mesh shader
             //context.startTimingQuery(context.getCommandList(PBMPM_G2P2G_COMPUTE_ID));
             scene->compute(renderModeType != 2);

        	 openXR->PollSystemEvents();
             openXR->PollEvents();

             if (openXR->IsSessionRunning()) {
                 // Render Frame
                 //context.startTimingQuery(context.getCommandList(PBMPM_G2P2G_COMPUTE_ID));
                 openXR->RenderFrame(*scene);
                 if (pbmpmIterConstants.mouseActivation == 1 || !PBMPMScene::constantsEqual(pbmpmIterConstants, pbmpmCurrConstants)) {
                     scene->updatePBMPMConstants(pbmpmIterConstants);
                     pbmpmCurrConstants = pbmpmIterConstants;
                 }
             }

                 //get pipelines
        /*
        auto renderPipeline = scene->getObjectSolidPipeline();
        //      auto fluidMeshPipeline = scene.getFluidMeshPipeline();
              //auto elasticMeshPipeline = scene.getElasticMeshPipeline();
              //auto viscoMeshPipeline = scene.getViscoMeshPipeline();
              //auto sandMeshPipeline = scene.getSandMeshPipeline();
              ////auto snowMeshPipeline = scene.getSnowMeshPipeline();
        //      auto objectWirePipeline = scene.getObjectWirePipeline();
        //      auto objectSolidPipeline = scene.getObjectSolidPipeline();
              //whichever pipeline renders first should begin and end the frame
              //auto firstPipeline = objectWirePipeline;

              //begin frame
        Window::get().beginFrame(renderPipeline->getCommandList());

        //create viewport
        D3D12_VIEWPORT vp;
        Window::get().createViewport(vp, renderPipeline->getCommandList());

        //wire object render pass
        Window::get().setRT(renderPipeline->getCommandList());
        Window::get().setViewport(vp, renderPipeline->getCommandList());
        //if (renderGrid) scene.drawGrid();
        //if (renderSpawn) scene.drawSpawners();
        scene->drawSolidObjects();

        //particles + imgui render pass
        //Window::get().setRT(renderPipeline->getCommandList());
        //Window::get().setViewport(vp, renderPipeline->getCommandList());
        // Only draw particles if we are not in the mesh shading mode
        if (renderModeType != 0) {
            scene->drawPBMPM();
        }

        //fluid mesh render pass
        if (scene->renderToggles[0]) {
            if (renderModeType != 2) scene->drawFluid(meshletRenderType, toonShadingLevels);
        }
        context->executeCommandList(renderPipeline->getCommandListID());
		context->resetCommandList(renderPipeline->getCommandListID());*/

        /*
        // elastic mesh render pass
        if (scene.renderToggles[1]) {
            Window::get().setRT(elasticMeshPipeline->getCommandList());
            Window::get().setViewport(vp, elasticMeshPipeline->getCommandList());
            if (renderModeType != 2) scene.drawElastic(meshletRenderType, toonShadingLevels);
            context.executeCommandList(elasticMeshPipeline->getCommandListID());
        }

        // sand mesh render pass
        if (scene.renderToggles[2]) {
            Window::get().setRT(sandMeshPipeline->getCommandList());
            Window::get().setViewport(vp, sandMeshPipeline->getCommandList());
            if (renderModeType != 2) scene.drawSand(meshletRenderType, toonShadingLevels);
            context.executeCommandList(sandMeshPipeline->getCommandListID());
        }

        // visco mesh render pass
        if (scene.renderToggles[3]) {
            Window::get().setRT(viscoMeshPipeline->getCommandList());
            Window::get().setViewport(vp, viscoMeshPipeline->getCommandList());
            if (renderModeType != 2) scene.drawVisco(meshletRenderType, toonShadingLevels);
            context.executeCommandList(viscoMeshPipeline->getCommandListID());
        }*/

        // snow mesh render pass
        /*if (scene.renderToggles[4]) {
            Window::get().setRT(snowMeshPipeline->getCommandList());
            Window::get().setViewport(vp, snowMeshPipeline->getCommandList());
            if (renderModeType != 2) scene.drawSnow(meshletRenderType, toonShadingLevels);
            context.executeCommandList(snowMeshPipeline->getCommandListID());
        }*/
        /*
        //set up ImGUI for frame
        Window::get().setRT(renderPipeline->getCommandList());
        Window::get().setViewport(vp, renderPipeline->getCommandList());
        ImGui_ImplDX12_NewFrame();
        ImGui_ImplWin32_NewFrame();
        ImGui::NewFrame();

        //draw ImGUI
        drawImGUIWindow(pbmpmIterConstants, io,
            scene->getFluidIsovalue(), 
            scene->getFluidKernelScale(), 
            scene->getFluidKernelRadius(),
        	scene->getElasticIsovalue(),
        	scene->getElasticKernelScale(),
        	scene->getElasticKernelRadius(),
        	scene->getSandIsovalue(),
        	scene->getSandKernelScale(),
        	scene->getSandKernelRadius(),
            scene->getViscoIsovalue(),
            scene->getViscoKernelScale(),
            scene->getViscoKernelRadius(),
            scene->getPBMPMSubstepCount(),
            scene->getNumParticles());

        //render ImGUI
        ImGui::Render();

        renderPipeline->getCommandList()->SetDescriptorHeaps(1, &imguiSRVHeap);
        ImGui_ImplDX12_RenderDrawData(ImGui::GetDrawData(), renderPipeline->getCommandList());

        context->executeCommandList(renderPipeline->getCommandListID());

        // reset the first pipeline so it can end the frame
        context->resetCommandList(renderPipeline->getCommandListID());
        //end frame
        Window::get().endFrame(renderPipeline->getCommandList());
        // Execute command list
        context->executeCommandList(renderPipeline->getCommandListID());

        Window::get().present();
        context->resetCommandList(renderPipeline->getCommandListID());*/
        //if (scene.renderToggles[0]) {
        //	context.resetCommandList(fluidMeshPipeline->getCommandListID());
        //}
  //      if (scene.renderToggles[1]) {
  //          context.resetCommandList(elasticMeshPipeline->getCommandListID());
  //      }
  //      if (scene.renderToggles[2]) {
        //	context.resetCommandList(sandMeshPipeline->getCommandListID());
        //}
        //if (scene.renderToggles[3]) {
        //	context.resetCommandList(viscoMeshPipeline->getCommandListID());
        //}
        ///*if (scene.renderToggles[4]) {
        //	context.resetCommandList(snowMeshPipeline->getCommandListID());
        //}*/

  //      context.resetCommandList(objectWirePipeline->getCommandListID());
  //      context.resetCommandList(objectSolidPipeline->getCommandListID());*/
    }

    // Now you can safely call ReportLiveDeviceObjects
    ComPointer<ID3D12DebugDevice> debugDevice;
    if (SUCCEEDED(context->getDevice()->QueryInterface(IID_PPV_ARGS(&debugDevice)))) {
    }

	ShutdownManager::Shutdown();

    if (debugDevice) {
        debugDevice->ReportLiveDeviceObjects(D3D12_RLDO_DETAIL);
    }

    /*ImGui_ImplDX12_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImGui::DestroyContext();
    imguiSRVHeap->Release();*/

    return 0;
}
