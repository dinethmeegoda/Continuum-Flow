#include "g2p2gRootSignature.hlsl"  // Includes the ROOTSIG definition
#include "../../pbmpmShaders/PBMPMCommon.hlsl"  // Includes the TileDataSize definition

// Taken from https://github.com/electronicarts/pbmpm

// Root constants bound to b0
ConstantBuffer<PBMPMConstants> g_simConstants : register(b0);

// Structured Buffer for particles (read-write UAV)
RWStructuredBuffer<Particle> g_particles : register(u0);

// Structured Buffer for free indices with atomic access (read-write UAV)
RWStructuredBuffer<int> g_freeIndices : register(u1);

// Structured Buffer for bukkit particle indices (read-only SRV)
StructuredBuffer<uint> g_bukkitParticleData : register(t0);

// Structured Buffer for bukkit thread data (read-only SRV)
StructuredBuffer<BukkitThreadData> g_bukkitThreadData : register(t1);

// Structured Buffer for grid source data (read-only SRV)
StructuredBuffer<int> g_gridSrc : register(t2);

// Structured Buffer for grid destination data (read-write UAV with atomic support)
RWStructuredBuffer<int> g_gridDst : register(u2);

// Structured Buffer for grid cells to be cleared (read-write UAV)
RWStructuredBuffer<int> g_gridToBeCleared : register(u3);

RWStructuredBuffer<int> g_tempTileData : register(u4);

RWStructuredBuffer<float4> g_positions : register(u5);

//groupshared int s_tileData[TileDataSize];
groupshared int s_tileDataDst[TileDataSize];

unsigned int localGridIndex(uint3 index) {
	return (index.z * TotalBukkitEdgeLength * TotalBukkitEdgeLength + index.y * TotalBukkitEdgeLength + index.x) * 5;
}

// Function to clamp a particle's position inside the guardian region of the grid
float3 projectInsideGuardian(float3 p, uint3 gridSize, float guardianSize)
{
    // Define the minimum and maximum clamp boundaries
    float3 clampMin = float3(guardianSize, guardianSize, guardianSize);
    float3 clampMax = float3(gridSize) - float3(guardianSize, guardianSize, guardianSize) - float3(1.0, 1.0, 1.0);

    // Clamp the position `p` to be within the defined boundaries
    return clamp(p, clampMin, clampMax);
}

// Matrix Helper Functions

// Structure to hold the SVD result

// Function to compute the determinant of a  3x3 matrix
float det(float3x3 m) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

// Function to compute the trace of a 2x2 matrix
float tr(float2x2 m)
{
    return m[0][0] + m[1][1];
}

// Function to compute the trace of a 3x3 matrix
float tr3D(float3x3 m) {
    return m[0][0] + m[1][1] + m[2][2];
}

// Function to compute the inverse of a 3x3 matrix
float3x3 inverse(float3x3 m) {
    
    float d = det(m);
    if (abs(d) < 1e-12) {
        return Identity;
    }

    float3x3 adj;
    adj[0][0] = +(m[1][1] * m[2][2] - m[2][1] * m[1][2]);
    adj[0][1] = -(m[0][1] * m[2][2] - m[2][1] * m[0][2]);
    adj[0][2] = +(m[0][1] * m[1][2] - m[1][1] * m[0][2]);
    adj[1][0] = -(m[1][0] * m[2][2] - m[2][0] * m[1][2]);
    adj[1][1] = +(m[0][0] * m[2][2] - m[2][0] * m[0][2]);
    adj[1][2] = -(m[0][0] * m[1][2] - m[1][0] * m[0][2]);
    adj[2][0] = +(m[1][0] * m[2][1] - m[2][0] * m[1][1]);
    adj[2][1] = -(m[0][0] * m[2][1] - m[2][0] * m[0][1]);
    adj[2][2] = +(m[0][0] * m[1][1] - m[1][0] * m[0][1]);
    return adj * (1.0 / d);
}

float3x3 outerProduct(float3 x, float3 y) {
    return float3x3(
        x.x * y.x, x.x * y.y, x.x * y.z,
        x.y * y.x, x.y * y.y, x.y * y.z,
        x.z * y.x, x.z * y.y, x.z * y.z
    );
}

// Function to create a diagonal 3x3 matrix from a float3 vector
float3x3 diag(float3 d)
{
    return float3x3(
        d.x, 0, 0,
        0, d.y, 0,
        0, 0, d.z
    );
}

// Function to truncate 4x4 matrix to 2x2 matrix
float2x2 truncate(float4x4 m)
{
    return float2x2(m[0].xy, m[1].xy);
}

struct SVDResult
{
    float3x3 U;
    float3 Sigma;
    float3x3 Vt;
};

float3x3 getRotationMatrix(float c, float s, int i, int j) {
    float3x3 R = Identity;
    R[i][i] = c;
    R[i][j] = -s;
    R[j][i] = s;
    R[j][j] = c;
    return R;
}

// Compute off-diagonal sum of squares
float offDiagonalSum(float3x3 mat) {
    float sum = 0.0f;
    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            sum += mat[i][j] * mat[i][j] + mat[j][i] * mat[j][i];
        }
    }
    return sum;
}

SVDResult svd(float3x3 A) {
    SVDResult result;
    const int MAX_ITERATIONS = 20; // Slightly more iterations
    const float EPSILON = 1e-6f;  // Slightly looser epsilon

    float3x3 U = Identity;
    float3x3 V = Identity;
    float3x3 B = A;

    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        bool converged = true;

        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                float3 col_i = float3(B[0][i], B[1][i], B[2][i]);
                float3 col_j = float3(B[0][j], B[1][j], B[2][j]);

                float a = dot(col_i, col_i);
                float c = dot(col_j, col_j);
                float b = dot(col_i, col_j);

                if (abs(b) < EPSILON * sqrt(a * c))
                    continue;

                converged = false;
                float zeta = (c - a) / (2.0f * b);
                float t = sign(zeta) / (abs(zeta) + sqrt(1.0f + zeta * zeta));
                float c_rot = 1.0f / sqrt(1.0f + t * t);
                float s_rot = c_rot * t;

                float3x3 J = getRotationMatrix(c_rot, s_rot, i, j);
                B = mul(B, J);
                V = mul(V, J);
            }
        }

        if (converged)
            break;
    }

    float3 singularValues;
    for (int i = 0; i < 3; i++) {
        float norm = length(float3(B[0][i], B[1][i], B[2][i]));
        singularValues[i] = norm > EPSILON ? norm : EPSILON; // Ensure no zero singular values

        if (norm > EPSILON) {
            B[0][i] /= norm;
            B[1][i] /= norm;
            B[2][i] /= norm;
        }
    }

    if (det(B) < 0) {
        B[0][2] = -B[0][2];
        B[1][2] = -B[1][2];
        B[2][2] = -B[2][2];
        singularValues[2] = -singularValues[2];
    }

    result.U = B;
    // Clamp singular values to prevent extreme distortion
    result.Sigma = clamp(singularValues, float3(0.5, 0.5, 0.5), float3(5000.0, 5000.0, 5000.0));
    result.Vt = transpose(V);

    return result;
}

[numthreads(ParticleDispatchSize, 1, 1)]
void main(uint indexInGroup : SV_GroupIndex, uint3 groupId : SV_GroupID)
{
    // Load thread-specific data
    BukkitThreadData threadData = g_bukkitThreadData[groupId.x];

    // Calculate grid origin
    int3 localGridOrigin = BukkitSize * int3(threadData.bukkitX, threadData.bukkitY, threadData.bukkitZ)
        - int3(BukkitHaloSize, BukkitHaloSize, BukkitHaloSize);
    int3 idInGroup = int3(
        indexInGroup % TotalBukkitEdgeLength,
        int(indexInGroup / TotalBukkitEdgeLength) % TotalBukkitEdgeLength,
        int(indexInGroup / (TotalBukkitEdgeLength * TotalBukkitEdgeLength)));
   
    int3 gridVertex = idInGroup + localGridOrigin;
    float3 gridPosition = float3(gridVertex);

    // Initialize variables
    float dx = 0.0;
    float dy = 0.0;
    float dz = 0.0;
    float w = 0.0; //weight
    float v = 0.0; //volume

    // Check if grid vertex is within valid bounds
    bool gridVertexIsValid = all(gridVertex >= int3(0, 0, 0)) && all(gridVertex <= g_simConstants.gridSize);

    if (gridVertexIsValid)
    {
        uint gridVertexAddress = gridVertexIndex(uint3(gridVertex), g_simConstants.gridSize);

		// Load grid data
        dx = decodeFixedPoint(g_gridSrc[gridVertexAddress + 0], g_simConstants.fixedPointMultiplier);
        dy = decodeFixedPoint(g_gridSrc[gridVertexAddress + 1], g_simConstants.fixedPointMultiplier);
        dz = decodeFixedPoint(g_gridSrc[gridVertexAddress + 2], g_simConstants.fixedPointMultiplier);
        w = decodeFixedPoint(g_gridSrc[gridVertexAddress + 3], g_simConstants.fixedPointMultiplier);
        v = decodeFixedPoint(g_gridSrc[gridVertexAddress + 4], g_simConstants.fixedPointMultiplier);

        // Grid update
        if (w < 1e-5f)
        {
            dx = 0.0f;
            dy = 0.0f;
            dz = 0.0f;
        }
        else
        {
            dx /= w;
            dy /= w;
            dz /= w;
        }

        float3 gridDisplacement = float3(dx, dy, dz);

        // Collision detection against guardian shape

        // Grid vertices near or inside the guardian region should have their displacement values
        // corrected in order to prevent particles moving into the guardian.
        // We do this by finding whether a grid vertex would be inside the guardian region after displacement
        // with the current velocity and, if it is, setting the displacement so that no further penetration can occur.

        float3 displacedGridPosition = gridPosition + gridDisplacement;
        float3 projectedGridPosition = projectInsideGuardian(displacedGridPosition, g_simConstants.gridSize, GuardianSize + 1);
        float3 projectedDifference = projectedGridPosition - displacedGridPosition;

        if (any(projectedDifference != 0))
        {
            // Calculate normal direction
            float3 normal = normalize(projectedDifference);
            // Project out the normal component
            float3 tangential = gridDisplacement - normal * dot(gridDisplacement, normal);
            // Apply friction to tangential component
            gridDisplacement = tangential * (1.0 - g_simConstants.borderFriction);
        }

        /*if (projectedDifference.x != 0)
        {
            gridDisplacement.x = 0;
            gridDisplacement.y = lerp(gridDisplacement.y, 0, g_simConstants.borderFriction);
            gridDisplacement.z = lerp(gridDisplacement.z, 0, g_simConstants.borderFriction);
        }

        if (projectedDifference.y != 0)
        {
            gridDisplacement.y = 0;
            gridDisplacement.x = lerp(gridDisplacement.x, 0, g_simConstants.borderFriction);
            gridDisplacement.z = lerp(gridDisplacement.z, 0, g_simConstants.borderFriction);
        }

        if (projectedDifference.z != 0)
        {
            gridDisplacement.z = 0;
            gridDisplacement.x = lerp(gridDisplacement.x, 0, g_simConstants.borderFriction);
            gridDisplacement.y = lerp(gridDisplacement.y, 0, g_simConstants.borderFriction);
        }*/

        dx = gridDisplacement.x;
        dy = gridDisplacement.y;
        dz = gridDisplacement.z;
    }

    // Save grid to local memory
    unsigned int tileDataIndex = localGridIndex(idInGroup);
    // Store encoded fixed-point values atomically
    int originalValue;
    //InterlockedExchange(s_tileData[tileDataIndex], encodeFixedPoint(dx, g_simConstants.fixedPointMultiplier), originalValue);
    //InterlockedExchange(s_tileData[tileDataIndex + 1], encodeFixedPoint(dy, g_simConstants.fixedPointMultiplier), originalValue);
    //InterlockedExchange(s_tileData[tileDataIndex + 2], encodeFixedPoint(dz, g_simConstants.fixedPointMultiplier), originalValue);
    //InterlockedExchange(s_tileData[tileDataIndex + 3], encodeFixedPoint(w, g_simConstants.fixedPointMultiplier), originalValue);
    //InterlockedExchange(s_tileData[tileDataIndex + 4], encodeFixedPoint(v, g_simConstants.fixedPointMultiplier), originalValue);
    InterlockedExchange(g_tempTileData[(groupId.x * TileDataSize) + tileDataIndex], encodeFixedPoint(dx, g_simConstants.fixedPointMultiplier), originalValue);
    InterlockedExchange(g_tempTileData[(groupId.x * TileDataSize) + tileDataIndex + 1], encodeFixedPoint(dy, g_simConstants.fixedPointMultiplier), originalValue);
    InterlockedExchange(g_tempTileData[(groupId.x * TileDataSize) + tileDataIndex + 2], encodeFixedPoint(dz, g_simConstants.fixedPointMultiplier), originalValue);
    InterlockedExchange(g_tempTileData[(groupId.x * TileDataSize) + tileDataIndex + 3], encodeFixedPoint(w, g_simConstants.fixedPointMultiplier), originalValue);
    InterlockedExchange(g_tempTileData[(groupId.x * TileDataSize) + tileDataIndex + 4], encodeFixedPoint(v, g_simConstants.fixedPointMultiplier), originalValue);
    
    // Make sure all values in destination grid are 0
    s_tileDataDst[tileDataIndex] = 0;
    s_tileDataDst[tileDataIndex + 1] = 0;
    s_tileDataDst[tileDataIndex + 2] = 0;
    s_tileDataDst[tileDataIndex + 3] = 0;
    s_tileDataDst[tileDataIndex + 4] = 0;

    // Synchronize all threads in the group
    GroupMemoryBarrierWithGroupSync();
    
    if (indexInGroup < threadData.rangeCount)
    {
        // Load Particle
        uint myParticleIndex = g_bukkitParticleData[threadData.rangeStart + indexInGroup];
        
        Particle particle = g_particles[myParticleIndex];
        float liquidDensity = g_positions[myParticleIndex].w;
        
        float3 p = g_positions[myParticleIndex].xyz;
        QuadraticWeightInfo weightInfo = quadraticWeightInit(p);
        
        if (g_simConstants.iteration != 0)
        {
            // G2P
            float3x3 B = ZeroMatrix;
            float3 d = float3(0, 0, 0);
            float volume = 0.0;
            
            // Iterate over local 3x3 neighborhood
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++) {
                        // Weight corresponding to this neighborhood cell
                        float weight = weightInfo.weights[i].x * weightInfo.weights[j].y * weightInfo.weights[k].z;

                        // Grid vertex index
                        int3 neighborCellIndex = int3(weightInfo.cellIndex) + int3(i, j, k);

                        // 3D index relative to the corner of the local grid
                        int3 neighborCellIndexLocal = neighborCellIndex - localGridOrigin;

                        // Linear Index in the local grid
                        uint gridVertexIdx = localGridIndex(uint3(neighborCellIndexLocal));

                        int fixedPoint0;
                        //InterlockedAdd(s_tileData[gridVertexIdx + 0], 0, fixedPoint0);
                        InterlockedAdd(g_tempTileData[(groupId.x * TileDataSize) + gridVertexIdx + 0], 0, fixedPoint0);
                        int fixedPoint1;
                        //InterlockedAdd(s_tileData[gridVertexIdx + 1], 0, fixedPoint1);
                        InterlockedAdd(g_tempTileData[(groupId.x * TileDataSize) + gridVertexIdx + 1], 0, fixedPoint1);
                        int fixedPoint2;
                        //InterlockedAdd(s_tileData[gridVertexIdx + 2], 0, fixedPoint2);
                        InterlockedAdd(g_tempTileData[(groupId.x * TileDataSize) + gridVertexIdx + 2], 0, fixedPoint2);

                        float3 weightedDisplacement = weight * float3(
                            decodeFixedPoint(fixedPoint0, g_simConstants.fixedPointMultiplier),
                            decodeFixedPoint(fixedPoint1, g_simConstants.fixedPointMultiplier),
                            decodeFixedPoint(fixedPoint2, g_simConstants.fixedPointMultiplier));

                        float3 offset = float3(neighborCellIndex) - p + 0.5;
                        B += outerProduct(weightedDisplacement, offset);
                        d += weightedDisplacement;

                        if (g_simConstants.useGridVolumeForLiquid != 0)
                        {
                            int fixedPoint4;
                            //InterlockedAdd(s_tileData[gridVertexIdx + 4], 0, fixedPoint4);
                            InterlockedAdd(g_tempTileData[(groupId.x * TileDataSize) + gridVertexIdx + 4], 0, fixedPoint4);
                            volume += weight * decodeFixedPoint(fixedPoint4, g_simConstants.fixedPointMultiplier);
                        }
                    }
                }

            }
            
            if (g_simConstants.useGridVolumeForLiquid != 0)
            {
                // Update particle volume
                
                float safeVolume = max(volume, 1e-6);
                volume = 1.0 / safeVolume;
                if (volume < 1.0)
                {
                    liquidDensity = lerp(liquidDensity, volume, 0.1);
                }
            }
            
            // Save the deformation gradient as a 3x3 matrix by adding the identity matrix to the rest
            particle.deformationDisplacement = B * 4.0;
            particle.displacement = d;
            
            // Integration
            if (g_simConstants.iteration == g_simConstants.iterationCount - 1)
            {
                if (particle.material == MaterialLiquid)
                {
                    // The liquid material only cares about the determinant of the deformation gradient.
                    // We can use the regular MPM integration below to evolve the deformation gradient, but
                    // this approximation actually conserves more volume.
                    // This is based on det(F^n+1) = det((I + D)*F^n) = det(I+D)*det(F^n)
                    // and noticing that D is a small velocity, we can use the identity det(I + D) ≈ 1 + tr(D) to first order
                    // ending up with det(F^n+1) = (1+tr(D))*det(F^n)
                    // Then we directly set particle.liquidDensity to reflect the approximately integrated volume.
                    // The liquid material does not actually use the deformation gradient matrix.
                    liquidDensity *= (tr3D(particle.deformationDisplacement) + 1.0);

                    // Safety clamp to avoid instability with very small densities.
                    liquidDensity = max(liquidDensity, 0.05);
                }
                if (particle.material != MaterialLiquid) {


                    SVDResult svdResult = svd(particle.deformationGradient);
                    // Safety clamp to prevent numerical instability
                    // Clamp each singular value to prevent extreme deformation
                    // 
                    svdResult.Sigma = clamp(svdResult.Sigma, float3(0.1, 0.1, 0.1), float3(10000.0, 10000.0, 10000.0));

                    if (particle.material == MaterialSand) {
                        // Drucker - Prager sand based on :
                        // Gergely Klár, Theodore Gast, Andre Pradhana, Chuyuan Fu, Craig Schroeder, Chenfanfu Jiang, and Joseph Teran. 2016.
                        // Drucker-prager elastoplasticity for sand animation. ACM Trans. Graph. 35, 4, Article 103 (July 2016), 12 pages.
                        // https://doi.org/10.1145/2897824.2925906
                        float sinPhi = sin(g_simConstants.frictionAngle * 3.14159 / 180.0);
                        float alpha = sqrt(2.0 / 3.0) * 2.0 * sinPhi / (3.0 - sinPhi);
                        float beta = 0.5;

                        float3 eDiag = log(max(abs(svdResult.Sigma), float3(1e-6, 1e-6, 1e-6)));
                        float3x3 eps = diag(eDiag);
                        float trace = tr3D(eps) + particle.logJp;

                        float3x3 eHat = eps - (trace / 3.0) * Identity;  // Note: Changed from 2 to 3 for 3D
                        float frobNrm = sqrt(dot(eHat[0], eHat[0]) +
                            dot(eHat[1], eHat[1]) +
                            dot(eHat[2], eHat[2]));



                        float elasticityRatio = 0.9f;
                        if (trace >= 0.0) {
                            // Instead of fully resetting, only partially relax towards identity
                            float relaxFactor = 0.5;
                            svdResult.Sigma = lerp(svdResult.Sigma, float3(1.0, 1.0, 1.0), relaxFactor);
                            particle.logJp = lerp(particle.logJp, beta * trace, 0.5);
                        }
                        else {
                            particle.logJp = 0;
                            float deltaGammaI = frobNrm + (elasticityRatio + 1.0) * trace * alpha;
                            if (deltaGammaI > 0) {
                                // Project to yield surface with a bit more damping
                                float damp = 0.5;
                                float3 h = eDiag - damp * deltaGammaI / frobNrm *
                                    (eDiag - float3(trace / 3.0, trace / 3.0, trace / 3.0));
                                svdResult.Sigma = exp(h);
                            }
                        }
                       
                    }
                    particle.deformationGradient = expandToFloat4x4(mul(mul(svdResult.U, diag(svdResult.Sigma)), svdResult.Vt));
                }
                
                // Update particle position
                p += particle.displacement;
                
                // Mouse Iteraction Here
                /*if (g_simConstants.mouseActivation > 0)
                {
                    float3 offset = particle.position - g_simConstants.mousePosition;
                    float lenOffset = max(length(offset), 0.0001);
                    if (lenOffset < g_simConstants.mouseRadius)
                    {
                        float3 normOffset = offset / lenOffset;

                        if (g_simConstants.mouseFunction == 0)
                        {
                            particle.displacement += normOffset * 500.f;
                        }
                        else if (g_simConstants.mouseFunction == 1)
                        {
                            particle.displacement = g_simConstants.mouseVelocity * g_simConstants.deltaTime;
                        }
                    }
                }*/
                
                // Gravity Acceleration is normalized to the vertical size of the window
                particle.displacement.y -= float(g_simConstants.gridSize.y) * g_simConstants.gravityStrength * g_simConstants.deltaTime * g_simConstants.deltaTime;
                
                // Free count may be negative because of emission. So make sure it is at last zero before incrementing.
                int originalMax; // Needed for InterlockedMax output parameter
                InterlockedMax(g_freeIndices[0], 0, originalMax); 
                
                p = projectInsideGuardian(p, g_simConstants.gridSize, GuardianSize);
            }
            
            // Save the particle back to the buffer
            g_particles[myParticleIndex] = particle;
			g_positions[myParticleIndex] = float4(p, liquidDensity);
        }
        
        {
            // Particle update
            if (particle.material == MaterialLiquid)
            {
                // Simple liquid viscosity: just remove deviatoric part of the deformation displacement
                float3x3 deviatoric = -1.0 * (particle.deformationDisplacement + transpose(particle.deformationDisplacement));
                particle.deformationDisplacement += g_simConstants.liquidViscosity * 0.5 * deviatoric;

                float alpha = 0.5 * (1.0 / liquidDensity - tr3D(particle.deformationDisplacement) - 1.0);
                particle.deformationDisplacement += g_simConstants.liquidRelaxation * alpha * Identity;
            }
            else if (particle.material == MaterialSand)
            {
                // Compute updated deformation gradient
                float3x3 F = mul((Identity + particle.deformationDisplacement), particle.deformationGradient);

                SVDResult svdResult = svd(F);

                // Drucker-Prager parameters
                float sinPhi = sin(g_simConstants.frictionAngle * 3.14159f / 180.0f);
                float alpha = sqrt(2.0f / 3.0f) * 2.0f * sinPhi / (3.0f - sinPhi);
                float beta = 0.5f;

                // Ensure no zero singular values
                float3 safeSigma = max(svdResult.Sigma, float3(1e-6f, 1e-6f, 1e-6f));
                float3 eDiag = float3(log(safeSigma.x), log(safeSigma.y), log(safeSigma.z));

                float3x3 eps = diag(eDiag);

                float trace = tr3D(eps) + particle.logJp;

                // Deviatoric part: subtract (trace/3)*Identity
                float3x3 eHat = eps - (trace / 3.0f) * Identity;

                // Frobenius norm of eHat
                float frobNrm = 0.0f;
                [unroll]
                    for (int row = 0; row < 3; row++)
                    {
                        [unroll]
                            for (int col = 0; col < 3; col++)
                            {
                                frobNrm += eHat[row][col] * eHat[row][col];
                            }
                    }
                frobNrm = sqrt(frobNrm);

                float elasticityRatio = 0.9;
                float elasticRelaxation = 1.5;

                if (trace >= 0.0f)
                {
                    // Expansion: reset Sigma to identity
                    svdResult.Sigma = float3(1.0f, 1.0f, 1.0f);
                    particle.logJp = beta * trace;
                }
                else
                {
                    particle.logJp = 0;
                    float deltaGammaI = frobNrm + (elasticityRatio + 1.0f) * trace * alpha;
                    if (deltaGammaI > 0.0f)
                    {
                        // Project to yield surface
                        // Move eDiag towards mean to reduce deformation
                        float meanStrain = trace / 3.0f;
                        float3 meanVec = float3(meanStrain, meanStrain, meanStrain);
                        float3 eDiagDiff = eDiag - meanVec;
                        float3 h = eDiag - (deltaGammaI / frobNrm) * eDiagDiff;

                        svdResult.Sigma = float3(exp(h.x), exp(h.y), exp(h.z));
                    }
                }

                // Clamp Sigma to avoid extreme values
                svdResult.Sigma = clamp(svdResult.Sigma, float3(1.0f, 1.0f, 1.0f), float3(1000.0f, 1000.0f, 1000.0f));

                // Reconstruct deformation gradient after plastic correction
                float3x3 SigmaMat = diag(svdResult.Sigma);
                float3x3 newF = mul(mul(svdResult.U, SigmaMat), svdResult.Vt);

                // Volume preservation: enforce det == 1 by blending with Q
                float df = det(newF);
                float cdf = clamp(abs(df), 0.1f, 1.0f);
                float3x3 Q = (1.0f / (sign(df) * cbrt(cdf))) * F;

                // Blend between newF and Q
                float alphaVal = elasticityRatio;
                float3x3 tgt = alphaVal * newF + (1.0f - alphaVal) * Q;

                float3x3 invOldF = inverse(particle.deformationGradient);
                float3x3 diff = (mul(tgt, invOldF) - Identity) - particle.deformationDisplacement;
                particle.deformationDisplacement +=elasticRelaxation * diff;

                // Apply a deviatoric viscosity
                float3x3 dev = -1.0f * (particle.deformationDisplacement + transpose(particle.deformationDisplacement));
                particle.deformationDisplacement += g_simConstants.liquidViscosity * 0.5f * dev;
            }

            // P2G

            // Iterate over local 3x3 neighborhood
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++) 
                    {
                        // Weight corresponding to this neighborhood cell
                        float weight = weightInfo.weights[i].x * weightInfo.weights[j].y * weightInfo.weights[k].z;

                        // Grid vertex index
                        int3 neighborCellIndex = int3(weightInfo.cellIndex) + int3(i, j, k);

                        // 3D index relative to the corner of the local grid
                        int3 neighborCellIndexLocal = neighborCellIndex - localGridOrigin;

                        // Linear Index in the local grid
                        uint gridVertexIdx = localGridIndex(uint3(neighborCellIndexLocal));

                        // Update grid data
                        float3 offset = float3(neighborCellIndex) - p + 0.5;

                        float weightedMass = weight * particle.mass;
                        float3 momentum = weightedMass * (particle.displacement + mul(particle.deformationDisplacement, offset));

                        InterlockedAdd(s_tileDataDst[gridVertexIdx + 0], encodeFixedPoint(momentum.x, g_simConstants.fixedPointMultiplier));
                        InterlockedAdd(s_tileDataDst[gridVertexIdx + 1], encodeFixedPoint(momentum.y, g_simConstants.fixedPointMultiplier));
                        InterlockedAdd(s_tileDataDst[gridVertexIdx + 2], encodeFixedPoint(momentum.z, g_simConstants.fixedPointMultiplier));
                        InterlockedAdd(s_tileDataDst[gridVertexIdx + 3], encodeFixedPoint(weightedMass, g_simConstants.fixedPointMultiplier));

                        if (g_simConstants.useGridVolumeForLiquid != 0)
                        {
                            InterlockedAdd(s_tileDataDst[gridVertexIdx + 4], encodeFixedPoint(weight * particle.volume, g_simConstants.fixedPointMultiplier));
                        }
                    }
                }
            }
        }
    }
    
    // Synchronize all threads in the group
    GroupMemoryBarrierWithGroupSync();
    
    // Save Grid
    if (gridVertexIsValid)
    {
        uint gridVertexAddress = gridVertexIndex(uint3(gridVertex), g_simConstants.gridSize);

        // Atomic loads from shared memory using InterlockedAdd with 0
        int dxi = s_tileDataDst[tileDataIndex + 0];
        int dyi = s_tileDataDst[tileDataIndex + 1];
        int dzi = s_tileDataDst[tileDataIndex + 2];
        int wi = s_tileDataDst[tileDataIndex + 3];
        int vi = s_tileDataDst[tileDataIndex + 4];

    // Atomic adds to the destination buffer
        InterlockedAdd(g_gridDst[gridVertexAddress + 0], dxi);
        InterlockedAdd(g_gridDst[gridVertexAddress + 1], dyi);
        InterlockedAdd(g_gridDst[gridVertexAddress + 2], dzi);
        InterlockedAdd(g_gridDst[gridVertexAddress + 3], wi);
        InterlockedAdd(g_gridDst[gridVertexAddress + 4], vi);
    
    // Clear the entries in g_gridToBeCleared
        g_gridToBeCleared[gridVertexAddress + 0] = 0;
        g_gridToBeCleared[gridVertexAddress + 1] = 0;
        g_gridToBeCleared[gridVertexAddress + 2] = 0;
        g_gridToBeCleared[gridVertexAddress + 3] = 0;
        g_gridToBeCleared[gridVertexAddress + 4] = 0;

        g_tempTileData[(groupId.x * TileDataSize) + tileDataIndex + 0] = 0;
        g_tempTileData[(groupId.x * TileDataSize) + tileDataIndex + 1] = 0;
        g_tempTileData[(groupId.x * TileDataSize) + tileDataIndex + 2] = 0;
        g_tempTileData[(groupId.x * TileDataSize) + tileDataIndex + 3] = 0;
        g_tempTileData[(groupId.x * TileDataSize) + tileDataIndex + 4] = 0;
    }

}