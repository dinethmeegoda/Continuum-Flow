#define ROOTSIG \
"RootFlags(ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT)," \
"RootConstants(num32BitConstants=51, b0),"   /* For matrices and color */ \
"DescriptorTable(CBV(b1))"                   /* TODO: Change to RootCBV instead of Descriptor Table for model matrices buffer */