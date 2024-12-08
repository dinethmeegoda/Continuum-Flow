#include "PBMPMRootSignature.hlsl"

// Particle materials as an SRV at register t1
StructuredBuffer<int> materials : register(t1);

struct PSInput
{
    float4 Position : SV_Position; // Position passed from vertex shader
    uint InstanceID : INSTANCE_ID; // Instance ID passed from vertex shader
};

[RootSignature(ROOTSIG)]
float4 main(PSInput input) : SV_Target
{
    // Use the instance ID to retrieve the material
    int materialIndex = materials[input.InstanceID];

// Determine the material color based on the material index (example logic)
if (materialIndex == 0)
    return float4(1.0f, 0.0f, 0.0f, 1.0f); // Red
else if (materialIndex == 1)
    return float4(0.0f, 0.75f, 0.0f, 1.0f); // Green
else if (materialIndex == 2)
    return float4(0.0f, 0.0f, 1.0f, 1.0f); // Blue
else
    return float4(1.0f, 1.0f, 1.0f, 1.0f); // Default: White
}
