#include "Pipeline.h"

Pipeline::Pipeline(std::string rootSignatureShaderName, DXContext& context, CommandListID cmdID,
	DescriptorHeap* dHeap)
	: rootSignatureShader(rootSignatureShaderName), descriptorHeap(dHeap), cmdID(cmdID),
	cmdList(context.getCommandList(cmdID))
{
	//context.resetCommandList(cmdID);
	context.getDevice()->CreateRootSignature(0, rootSignatureShader.getBuffer(), rootSignatureShader.getSize(), IID_PPV_ARGS(&rootSignature));
}

ComPointer<ID3D12RootSignature>& Pipeline::getRootSignature()
{
	return this->rootSignature;
}

DescriptorHeap* Pipeline::getDescriptorHeap()
{
	return descriptorHeap;
}

void Pipeline::releaseResources()
{
	rootSignature.Release();
}