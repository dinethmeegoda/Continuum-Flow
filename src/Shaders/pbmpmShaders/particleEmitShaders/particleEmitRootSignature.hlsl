#define ROOTSIG \
"RootFlags(ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT)," \
"RootConstants(num32BitConstants=22, b0)," \
"RootConstants(num32BitConstants=12, b1)," \
"CBV(b2)," /* For Sim Shapes */ \
"DescriptorTable(UAV(u0, numDescriptors=3)),"    /* Table for particleBuffer, freeIndicesBuffer, particleCountBuffer */ \
"DescriptorTable(SRV(t0, numDescriptors=1)), " /* Table for curr grid */ \
"DescriptorTable(UAV(u3, numDescriptors=3))," /* Table for g_positions & materials & displacement */ \
"DescriptorTable(UAV(u6, numDescriptors=1))" /* Table for mass and volume */