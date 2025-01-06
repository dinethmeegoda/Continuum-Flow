#include "PBMPMRootSignature.hlsl"
// Particle materials as an SRV at register t1
StructuredBuffer<float4> materials : register(t1);
struct PSInput
{
    float4 Position : SV_Position; // Position passed from vertex shader
    uint InstanceID : INSTANCE_ID; // Instance ID passed from vertex shader
};

float2 rand_2(in float2 uv) {
    float noiseX = (frac(sin(dot(uv, float2(12.9898, 78.233) * 2.0)) * 43758.5453));
    float noiseY = sqrt(1 - noiseX * noiseX);
    return float2(noiseX, noiseY);
}

[RootSignature(ROOTSIG)]
float4 main(PSInput input) : SV_Target
{
    // Determine the material color based on the data from material buffer
    return float4(materials[input.InstanceID].xyz, 1.0f);
}