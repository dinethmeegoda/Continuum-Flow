#define ROOTSIG \
"RootFlags(ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT)," \
"RootConstants(num32BitConstants=20, b0),"     /* 17 constants: see PBMPM Constants */ \
"DescriptorTable(SRV(t0, numDescriptors=3)),"     /* SRV table for g_gridSrc, g_bukkitThreadData, g_bukkitParticleData */ \
"DescriptorTable(UAV(u0, numDescriptors=4))"      /* UAV table for g_particles, g_gridDst, g_gridToBeCleared, g_freeIndices */

