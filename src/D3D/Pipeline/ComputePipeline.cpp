#include "ComputePipeline.h"

ComputePipeline::ComputePipeline(std::string rootSignatureShaderName, const std::string shaderFilePath, DXContext& context,
	CommandListID cmdID, DescriptorHeap* dHeap)
	: Pipeline(rootSignatureShaderName, context, cmdID, dHeap),
	computeShader(shaderFilePath)
{
	createPSOD();
	createPipelineState(context.getDevice());
}

void ComputePipeline::createPSOD()
{
	psoDesc.pRootSignature = rootSignature.Get();
	psoDesc.CS = CD3DX12_SHADER_BYTECODE(computeShader.getBuffer(), computeShader.getSize());

}

void ComputePipeline::createPipelineState(ComPointer<ID3D12Device6>& device)
{
	HRESULT hr = device->CreateComputePipelineState(&psoDesc, IID_PPV_ARGS(&pso));
	if (FAILED(hr)) {
		throw std::runtime_error("Failed to create compute pipeline state");
	}
}