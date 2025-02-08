#include "particleEmitRootSignature.hlsl"  // Includes the ROOTSIG definition
#include "../../pbmpmShaders/PBMPMCommon.hlsl"  // Includes the TileDataSize definition

// Taken from https://github.com/electronicarts/pbmpm

// Root constants bound to b0 & b1
cbuffer simConstants : register(b0) {
	PBMPMConstants g_simConstants;
};

cbuffer mouseConstants : register(b1) {
    MouseConstants g_mouseConstants;
};

// Define the constant buffer with an array of SimShapes
cbuffer shapes : register(b2)
{
    SimShape g_shapes[MaxSimShapes]; // Adjust the size of the array as needed
};

// Structured Buffer for particles (read-write UAV)
RWStructuredBuffer<Particle> g_particles : register(u0);

// Structured Buffer for free indices with atomic access (read-write UAV)
RWStructuredBuffer<int> g_freeIndices : register(u1);

// Structured Buffer for free indices with atomic access (read-write UAV)
RWStructuredBuffer<int> g_particleCount : register(u2);

// Structured Buffer for grid source data (read-only SRV)
StructuredBuffer<int> g_grid : register(t0);

// Structured Buffer for positions (read-write UAV)
RWStructuredBuffer<float4> g_positions : register(u3);

// Structured Buffer for materials (read-write UAV), color in first three components, material enum stored in fourth
RWStructuredBuffer<float4> g_materials : register(u4);

// Structured Buffer for displacement (read-write UAV)
RWStructuredBuffer<float4> g_displacement : register(u5);

// Structured Buffer for mass and volume (read-write UAV)
RWStructuredBuffer<float4> g_massVolume : register(u6);

static const float3 colorTable[] = {
    float3(0.0, 0.573, 0.878), // Water
    float3(0.0, 0.75, 0.0), // Elastic
    float3(0.8, 0.8, 0.0), // Sand
    float3(0.7, 0.0, 0.8), // Visco
    float3(0.8, 0.8, 0.8), // Snow
    float3(0.0, 0.0, 0.0)  // Default
};

uint hash(uint input)
{
    uint state = input * 747796405 + 2891336453;
    uint word = ((state >> ((state >> 28) + 4)) ^ state) * 277803737;
    return (word >> 22) ^ word;
}

float3 generateJitter(float3 seed)
{
    // Scale factor to map the range [0, 1] to [-0.25, 0.25]
    const float scale = 0.25;

    // Use a simple hash-like function to generate values from the seed
    float3 hashed = frac(sin(dot(seed, float3(12.9898, 78.233, 37.719))) * 43758.5453);

    // Map the range [0, 1] to [-0.25, 0.25]
    return (hashed * 2.0 - 1.0) * scale;
}


bool insideGuardian(uint3 id, uint3 gridSize, uint guardianSize)
{
    if (id.x <= guardianSize) { return false; }
    if (id.x >= (gridSize.x - guardianSize - 1)) { return false; }
    if (id.y <= guardianSize) { return false; }
    if (id.y >= gridSize.y - guardianSize - 1) { return false; }
    if (id.z <= guardianSize) { return false; }
    if (id.z >= gridSize.z - guardianSize - 1) { return false; }

    return true;
}

Particle createParticle()
{
    Particle particle;

    particle.deformationGradient = Identity;
    particle.deformationDisplacement = ZeroMatrix;

    particle.lambda = 0.0;
    particle.logJp = 1.0;
    particle.enabled = 1.0;

    return particle;
}

void addParticle(float3 position, int material, float volume, float density, float jitterScale)
{
	// If we're over the particle limit, don't add any more
    uint particleIndex = 0;
    // First check the free list to see if we can reuse a particle slot
    int freeIndexSlot;
    InterlockedAdd(g_freeIndices[0], -1, freeIndexSlot);
    freeIndexSlot--;

    if (freeIndexSlot >= 0)
    {
        InterlockedAdd(g_freeIndices[freeIndexSlot + 1], 0, particleIndex);
    }
    else // If free list is empty then grow the particle count
    {
        InterlockedAdd(g_particleCount[0], 1, particleIndex);
    }

    /*uint jitterX = hash(particleIndex);
    uint jitterY = hash(uint(position.x * position.y * 100));
    uint jitterZ = hash(uint(position.x * position.z * 200));*/

	float3 jitter = generateJitter(position);

    Particle newParticle = createParticle();

    g_particles[particleIndex] = newParticle;
    g_materials[particleIndex] = float4(colorTable[material], material);
	g_positions[particleIndex] = float4(position + jitter * jitterScale, 1.0);
	g_displacement[particleIndex] = float4(0, 0, 0, 0);
	g_massVolume[particleIndex] = float4(volume * density, volume, 0, 0);
}

[numthreads(GridDispatchSize, GridDispatchSize, GridDispatchSize)]
void main(uint3 id : SV_DispatchThreadID)
{
    if (!insideGuardian(id.xyz, g_simConstants.gridSize, GuardianSize + 1) || g_particleCount[0] >= MaxParticles)
    {
        return;
    }

    float3 pos = float3(id.xyz);

    QuadraticWeightInfo weightInfo = quadraticWeightInit(pos);
    int3 nearestCell = int3(weightInfo.cellIndex) + int3(1, 1, 1);
    float nearestCellVolume = decodeFixedPoint(g_grid[gridVertexIndex(uint3(nearestCell), g_simConstants.gridSize) + 4], g_simConstants.fixedPointMultiplier);

    for (int shapeIndex = 0; shapeIndex < g_simConstants.shapeCount; shapeIndex++)
    {
        SimShape shape = g_shapes[shapeIndex];

        bool isEmitter = shape.functionality == ShapeFunctionEmit;
        bool isInitialEmitter = shape.functionality == ShapeFunctionInitialEmit;

        if (!(isEmitter || isInitialEmitter))
        {
            continue;
        }

        // Skip emission if we are spewing liquid into an already compressed space
        if (isEmitter && shape.material == MaterialLiquid && nearestCellVolume > 1.5)
        {
            continue;
        }

        uint particleCountPerCellAxis = uint(g_simConstants.particlesPerCellAxis);
        float volumePerParticle = 1.0f / float(particleCountPerCellAxis * particleCountPerCellAxis);

        CollideResult c = collide(shape, pos);
        if (c.collides)
        {
            uint emitEvery = uint(1.0 / (shape.emissionRate * g_simConstants.deltaTime));

            for (int i = 0; i < particleCountPerCellAxis; i++)
            {
                for (int j = 0; j < particleCountPerCellAxis; j++)
                {
                    for (int k = 0; k < particleCountPerCellAxis; k++)
                    {
                        uint hashCodeX = hash(id.x * particleCountPerCellAxis + i);
                        uint hashCodeY = hash(id.y * particleCountPerCellAxis + j);
						uint hashCodeZ = hash(id.z * particleCountPerCellAxis + k);
						uint hashCode = hash(hashCodeX + hashCodeY + hashCodeZ);

                        bool emitDueToMyTurnHappening = isEmitter && 0 == ((hashCode + g_simConstants.simFrame) % emitEvery);
                        bool emitDueToInitialEmission = isInitialEmitter && g_simConstants.simFrame == 0;

                        float3 emitPos = pos + float3(float(i), float(j), float(k)) / float(particleCountPerCellAxis);

                        if (emitDueToInitialEmission || emitDueToMyTurnHappening)
                        {
                            addParticle(emitPos, shape.material, volumePerParticle, 1.0, 1.0 / float(particleCountPerCellAxis));
                        }
                    }
                }
            }
        }
    }
}