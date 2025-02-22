#pragma once

#include "../../Support/WinInclude.h"
#include "../../Support/Shader.h"
#include "../../Support/ComPointer.h"
#include "../../Support/Window.h"
#include "../DescriptorHeap.h"

class Pipeline {
public:
	Pipeline() = delete;
	Pipeline(std::string rootSignatureShaderName, DXContext& context, CommandListID cmdID,
		D3D12_DESCRIPTOR_HEAP_TYPE type, unsigned int numberOfDescriptors, D3D12_DESCRIPTOR_HEAP_FLAGS flags);
	~Pipeline() = default;

	virtual void createPSOD() = 0;
	virtual void createPipelineState(ComPointer<ID3D12Device6> device) = 0;

	ComPointer<ID3D12RootSignature>& getRootSignature();
	DescriptorHeap* getDescriptorHeap();
	ComPointer<ID3D12PipelineState>& getPSO() { return pso; }
	ID3D12GraphicsCommandList6* getCommandList() { return cmdList; }
	CommandListID getCommandListID() { return cmdID; }

	void releaseResources();

protected:
	Shader rootSignatureShader;

	ComPointer<ID3D12RootSignature> rootSignature;
	DescriptorHeap descriptorHeap;
	ComPointer<ID3D12PipelineState> pso;
	CommandListID cmdID;

	ID3D12GraphicsCommandList6* cmdList;
};