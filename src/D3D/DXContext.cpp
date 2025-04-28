#include "DXContext.h"
#include <iostream>

DXContext::DXContext() {

    if (FAILED(CreateDXGIFactory2(0, IID_PPV_ARGS(&dxgiFactory)))) {
        //handle could not create dxgi factory
        throw std::runtime_error("Could not create dxgi factory");
    }

    if (FAILED(D3D12CreateDevice(nullptr, D3D_FEATURE_LEVEL_11_0, IID_PPV_ARGS(&device)))) {
        //handle could not create device
        throw std::runtime_error("Could not create device");
    }

    D3D12_COMMAND_QUEUE_DESC cmdQueueDesc{};
    cmdQueueDesc.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;
    cmdQueueDesc.Priority = D3D12_COMMAND_QUEUE_PRIORITY_HIGH;
    cmdQueueDesc.NodeMask = 0;
    cmdQueueDesc.Flags = D3D12_COMMAND_QUEUE_FLAG_NONE;
    if (FAILED(device->CreateCommandQueue(&cmdQueueDesc, IID_PPV_ARGS(&cmdQueue)))) {
        //handle could not create command queue
        throw std::runtime_error("Could not create command queue");
    }

    if (FAILED(device->CreateFence(fenceValue, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&fence)))) {
        //handle could not create fence
        throw std::runtime_error("Could not create fence");
    }

    fenceEvent = CreateEvent(nullptr, false, false, nullptr);
    if (!fenceEvent) {
        //handle could not create fence event
        throw std::runtime_error("Could not create fence event");
    }

	for (int i = 0; i < NUM_CMDLISTS; i++) {
		ComPointer<ID3D12CommandAllocator> cmdAllocator;
        if (FAILED(device->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&cmdAllocator)))) {
			//handle cannot create cmd allocator
			throw std::runtime_error("Could not create command allocator");
		}
		cmdAllocators[i] = cmdAllocator;
        createCommandList(CommandListID(i));
		resetCommandList(CommandListID(i));
	}

    initTimingResources();

}

DXContext::~DXContext() {
    // 1. Make sure GPU is idle first
    flush(1); // or flush(FRAME_COUNT) if you have multiple frames

    // 2. Reset command lists (optional, but safest)
    for (auto& cmdList : cmdLists) {
        if (cmdList) {
            cmdList->Close();  // Close if not already closed (D3D12 expects closed command lists before release)
        }
    }

    // 3. Release command lists
    for (auto& cmdList : cmdLists) {
        if (cmdList) {
            cmdList.Release();
        }
    }

    // 4. Release command allocators
    for (auto& cmdAllocator : cmdAllocators) {
        if (cmdAllocator) {
            cmdAllocator.Release();
        }
    }

    // 5. Release query heap and buffer
    if (queryHeap) queryHeap.Release();
    if (queryResultBuffer) queryResultBuffer.Release();

    // 6. Release fence
    if (fence) fence.Release();

    // 7. Close fence event
    if (fenceEvent) {
        CloseHandle(fenceEvent);
        fenceEvent = nullptr;
    }

    // 8. Release device, command queue, factory
    if (cmdQueue) cmdQueue.Release();
    if (device) device.Release();
    if (dxgiFactory) dxgiFactory.Release();
}

void DXContext::signalAndWait() {
    cmdQueue->Signal(fence, ++fenceValue);
    if (SUCCEEDED(fence->SetEventOnCompletion(fenceValue, fenceEvent))) {
        if (WaitForSingleObject(fenceEvent, 20000) != WAIT_OBJECT_0) {
            std::exit(-1);
        }
    } else {
        std::exit(-1);
    }
}

void DXContext::resetCommandList(CommandListID id)
{
    cmdAllocators[id]->Reset();
    
	cmdLists[id]->Reset(cmdAllocators[id], nullptr);
}

void DXContext::executeCommandList(CommandListID id) {
	if (SUCCEEDED(cmdLists[id]->Close())) {
		ID3D12CommandList* lists[] = { cmdLists[id] };
		cmdQueue->ExecuteCommandLists(1, lists);
		signalAndWait();
	}
}

void DXContext::flush(size_t count) {
    for (size_t i = 0; i < count; i++) {
        signalAndWait();
    }
}

void DXContext::signalAndWaitForFence(ComPointer<ID3D12Fence>& fence, UINT64& fenceValue) {
	cmdQueue->Signal(fence, fenceValue);
    if (fence->GetCompletedValue() < fenceValue) {
        HANDLE eventHandle = CreateEvent(nullptr, FALSE, FALSE, nullptr);
        if (eventHandle == nullptr) {
            throw std::runtime_error("Failed to create event handle.");
        }

        // Set the event to be triggered when the GPU reaches the fence value
        fence->SetEventOnCompletion(fenceValue, eventHandle);

        // Wait until the event is triggered, meaning the GPU has finished
        WaitForSingleObject(eventHandle, INFINITE);
        CloseHandle(eventHandle);
    }
    fenceValue++;
}

void DXContext::initTimingResources() {
    // Create a query heap for timestamp queries
    D3D12_QUERY_HEAP_DESC queryHeapDesc = {};
    queryHeapDesc.Type = D3D12_QUERY_HEAP_TYPE_TIMESTAMP;
    queryHeapDesc.Count = 2; // One for start and one for end timestamp
    queryHeapDesc.NodeMask = 0;

    device->CreateQueryHeap(&queryHeapDesc, IID_PPV_ARGS(&queryHeap));

    // Create a resource to store the query results
    D3D12_RESOURCE_DESC bufferDesc = CD3DX12_RESOURCE_DESC::Buffer(2 * sizeof(UINT64));
    CD3DX12_HEAP_PROPERTIES heapProps(D3D12_HEAP_TYPE_READBACK);
    device->CreateCommittedResource(
        &heapProps,
        D3D12_HEAP_FLAG_NONE,
        &bufferDesc,
        D3D12_RESOURCE_STATE_COPY_DEST,
        nullptr,
        IID_PPV_ARGS(&queryResultBuffer));

}

double DXContext::readTimingQueryData() {
    UINT64* queryData;
    UINT64 startTime, endTime;
    queryResultBuffer->Map(0, nullptr, reinterpret_cast<void**>(&queryData));
    startTime = queryData[0];
    endTime = queryData[1];
    queryResultBuffer->Unmap(0, nullptr);

    UINT64 gpuFrequency;
    cmdQueue->GetTimestampFrequency(&gpuFrequency);
    double timeInMs = (endTime - startTime) / static_cast<double>(gpuFrequency) * 1000.0;

    return timeInMs;
}

void DXContext::startTimingQuery(ID3D12GraphicsCommandList6* cmdList) {
    cmdList->EndQuery(queryHeap.Get(), D3D12_QUERY_TYPE_TIMESTAMP, 0);
}

void DXContext::endTimingQuery(ID3D12GraphicsCommandList6* cmdList) {
    cmdList->EndQuery(queryHeap.Get(), D3D12_QUERY_TYPE_TIMESTAMP, 1);
    cmdList->ResolveQueryData(queryHeap.Get(), D3D12_QUERY_TYPE_TIMESTAMP, 0, 2, queryResultBuffer.Get(), 0);
}

ComPointer<IDXGIFactory7>& DXContext::getFactory() {
    return dxgiFactory;
}

ComPointer<ID3D12Device6>& DXContext::getDevice() {
    return device;
}

ComPointer<ID3D12CommandQueue>& DXContext::getCommandQueue() {
    return cmdQueue;
}

ID3D12GraphicsCommandList6* DXContext::createCommandList(CommandListID id)
{
    ComPointer<ID3D12GraphicsCommandList6> cmdList;
    if (FAILED(device->CreateCommandList1(0, D3D12_COMMAND_LIST_TYPE_DIRECT, D3D12_COMMAND_LIST_FLAG_NONE, IID_PPV_ARGS(&cmdList)))) {
        //handle could not create cmd list
        throw std::runtime_error("Could not create command list");
    }

    cmdLists[id] = cmdList;
    return cmdList.Get();
}
