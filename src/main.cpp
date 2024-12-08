#include "main.h"

int main() {
    //set up DX, window, keyboard mouse
    DebugLayer debugLayer = DebugLayer();
    DXContext context = DXContext();
    std::unique_ptr<Camera> camera = std::make_unique<Camera>();
    std::unique_ptr<Keyboard> keyboard = std::make_unique<Keyboard>();
    std::unique_ptr<Mouse> mouse = std::make_unique<Mouse>();

    if (!Window::get().init(&context, SCREEN_WIDTH, SCREEN_HEIGHT)) {
        //handle could not initialize window
        std::cout << "could not initialize window\n";
        Window::get().shutdown();
        return false;
    }

    //initialize ImGUI
    ImGuiIO& io = initImGUI(context);

    //set mouse to use the window
    mouse->SetWindow(Window::get().getHWND());

	// Get the client area of the window
    RECT rect;
    GetClientRect(Window::get().getHWND(), &rect);
    float clientWidth = static_cast<float>(rect.right - rect.left);
    float clientHeight = static_cast<float>(rect.bottom - rect.top);

    //initialize scene
    Scene scene{PBMPM, camera.get(), &context};

	//initialize pbmpm constants
	PBMPMConstants pbmpmConstants = scene.getPBMPMConstants();
    PBMPMConstants pbmpmTempConstants = pbmpmConstants;

    while (!Window::get().getShouldClose()) {
        //update window
        Window::get().update();
        if (Window::get().getShouldResize()) {
            //flush pending buffer operations in swapchain
            context.flush(FRAME_COUNT);
            Window::get().resize();
            camera->updateAspect((float)Window::get().getWidth() / (float)Window::get().getHeight());
        }

        //check keyboard state
        auto kState = keyboard->GetState();
        if (kState.W) {
            camera->translate({ 0.f, 0.f, 1.0f });
        }
        if (kState.A) {
            camera->translate({ -1.0f, 0.f, 0.f });
        }
        if (kState.S) {
            camera->translate({ 0.f, 0.f, -1.0f });
        }
        if (kState.D) {
            camera->translate({ 1.0f, 0.f, 0.f });
        }
        if (kState.Space) {
            camera->translate({ 0.f, 1.0f, 0.f });
        }
        if (kState.LeftControl) {
            camera->translate({ 0.f, -1.0f, 0.f });
        }
        if (kState.D1) {
            scene.setRenderScene(Object);
        }
        if (kState.D2) {
            scene.setRenderScene(PBMPM);
        }
        if (kState.D3) {
            scene.setRenderScene(Fluid);
        }

        //check mouse state
        auto mState = mouse->GetState();

        mouse->SetMode(mState.leftButton ? Mouse::MODE_RELATIVE : Mouse::MODE_ABSOLUTE);

        if (mState.positionMode == Mouse::MODE_RELATIVE && kState.LeftShift) {
            camera->rotateOnX(-mState.y * 0.01f);
            camera->rotateOnY(mState.x * 0.01f);
            camera->rotate();
        }

        if (mState.rightButton) {
            //enable mouse force
			bool previousMouseActivation = pbmpmTempConstants.mouseActivation > 0;

            pbmpmTempConstants.mouseActivation = 1;

            POINT mousePos;
            GetCursorPos(&mousePos);

			ScreenToClient(Window::get().getHWND(), &mousePos);

            float ndcX = (2.0f * mousePos.x / clientWidth) - 1.0f;
            float ndcY = 1.0f - (2.0f * mousePos.y / clientHeight);

			XMFLOAT4 prevMousePos = pbmpmTempConstants.mousePosition;
            pbmpmTempConstants.mousePosition = GetMouseWorldPositionAtDepth(Window::get().getHWND(), ndcX, ndcY, camera->getProjMat(), camera->getViewMat(), 70);

			// print the previous mouse position and the current mouse position
			std::cout << "Previous Mouse Position: " << prevMousePos.x << ", " << prevMousePos.y << ", " << prevMousePos.z << std::endl;
			std::cout << "Current Mouse Position: " << pbmpmTempConstants.mousePosition.x << ", " << pbmpmTempConstants.mousePosition.y << ", " << pbmpmTempConstants.mousePosition.z << std::endl;

            pbmpmTempConstants.mouseFunction = 0;
            pbmpmTempConstants.mouseRadius = 100;

			// If the mouse was previously activated, update the velocity
			if (previousMouseActivation) {
				pbmpmTempConstants.mouseVelocity = XMVectorGetX(XMVector3Length(XMLoadFloat4(&pbmpmTempConstants.mousePosition) - XMLoadFloat4(&prevMousePos))) * pbmpmTempConstants.deltaTime;
            }
            else {
				pbmpmTempConstants.mouseVelocity = pbmpmTempConstants.deltaTime;
			}

            scene.updatePBMPMConstants(pbmpmTempConstants);
        }
        else {
            pbmpmTempConstants.mouseActivation = 0;
			scene.updatePBMPMConstants(pbmpmTempConstants);
        }

        //update camera
        camera->updateViewMat();

        //get pipelines
        auto renderPipeline = scene.getRenderPipeline();
        auto meshPipeline = scene.getMeshPipeline();
        //whichever pipeline renders first should begin and end the frame
        auto firstPipeline = meshPipeline;

        //compute pbmpm + mesh shader
        scene.compute();

        //begin frame
        Window::get().beginFrame(firstPipeline->getCommandList());

        //create viewport
        D3D12_VIEWPORT vp;
        Window::get().createViewport(vp, firstPipeline->getCommandList());

        //mesh render pass
        Window::get().setRT(meshPipeline->getCommandList());
        Window::get().setViewport(vp, meshPipeline->getCommandList());
        scene.drawFluid();
        context.executeCommandList(meshPipeline->getCommandListID());

        //first render pass
        Window::get().setRT(renderPipeline->getCommandList());
        Window::get().setViewport(vp, renderPipeline->getCommandList());
        scene.draw();

        //set up ImGUI for frame
        ImGui_ImplDX12_NewFrame();
        ImGui_ImplWin32_NewFrame();
        ImGui::NewFrame();

        //draw ImGUI
        drawImGUIWindow(pbmpmTempConstants, io);

        //render ImGUI
        ImGui::Render();
        if (!PBMPMScene::constantsEqual(pbmpmTempConstants, pbmpmConstants)) {
            scene.updatePBMPMConstants(pbmpmTempConstants);
            pbmpmConstants = pbmpmTempConstants;
        }

        renderPipeline->getCommandList()->SetDescriptorHeaps(1, &imguiSRVHeap);
        ImGui_ImplDX12_RenderDrawData(ImGui::GetDrawData(), renderPipeline->getCommandList());

        context.executeCommandList(renderPipeline->getCommandListID());

        //end frame
        Window::get().endFrame(firstPipeline->getCommandList());

        Window::get().present();
		context.resetCommandList(renderPipeline->getCommandListID());
        context.resetCommandList(meshPipeline->getCommandListID());
    }

    // Scene should release all resources, including their pipelines
    scene.releaseResources();

    ImGui_ImplDX12_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImGui::DestroyContext();

    imguiSRVHeap->Release();

    //flush pending buffer operations in swapchain
    context.flush(FRAME_COUNT);
    Window::get().shutdown();
}
