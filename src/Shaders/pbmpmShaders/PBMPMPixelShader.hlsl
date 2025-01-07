#include "PBMPMRootSignature.hlsl"

cbuffer CameraMatrices : register(b0) {
    float4x4 viewMatrix;        // 16 floats
    float4x4 projectionMatrix;  // 16 floats
    float4x4 modelMatrix;       // 16 floats
    unsigned int renderMode;    // 1 int
};

// Particle positions as an SRV at register t0
StructuredBuffer<float4> positions : register(t0);

// Particle materials as an SRV at register t1
StructuredBuffer<float4> materials : register(t1);

static const float falloffTable[] = {
    0.75, // Water
    0.33, // Elastic
    0.75, // Sand
    0.2, // Visco
    0.7, // Snow
    0.5 // Default
};

struct PSInput
{
    float4 Position : SV_Position; // Position passed from vertex shader
    uint InstanceID : INSTANCE_ID; // Instance ID passed from vertex shader
};

float getBias(float time, float bias)
{
    return (time / ((((1.0 / bias) - 2.0) * (1.0 - time)) + 1.0));
}

[RootSignature(ROOTSIG)]
float4 main(PSInput input) : SV_Target
{
    // Get camera position from view matrix
	float3 cameraPosition = float3(-viewMatrix[0][3], -viewMatrix[1][3], -viewMatrix[2][3]);
    
    // Determine the material color based on the data from material buffer
    float3 color = materials[input.InstanceID].xyz;

    // Make the color darker further away from the camera
    float distance = length(cameraPosition - positions[input.InstanceID].xyz);
	// Make the fade distance based on the distance from the center of the grid, hard coded to 32x32x32 grid, change later
	float distanceFromGridCenter = length(cameraPosition - float3(16, 0, 16));
	float distanceRatio = clamp(distance * distance/(distanceFromGridCenter * 110), 0.0, 1.0);
    float3 finalColor = lerp(color * 1.1, color * falloffTable[materials[input.InstanceID].w], distanceRatio);

    // Return the final color
    return float4(finalColor, 1.0f);
}