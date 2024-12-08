#include "RootSignature.hlsl"

cbuffer CameraMatrices : register(b0) {
    float4x4 viewMatrix;        // 16 floats
    float4x4 projectionMatrix;  // 16 floats
    float4x4 modelMatrix;
    float3 color;
};

[RootSignature(ROOTSIG)]
float4 main() : SV_Target
{
    return float4(color.x, color.y, color.z, 1.0f);
}