#include "ifftw.h"

extern void X(codelet_hc2cfdftv_2)(planner *);
extern void X(codelet_hc2cfdftv_4)(planner *);
extern void X(codelet_hc2cfdftv_6)(planner *);
extern void X(codelet_hc2cfdftv_8)(planner *);
extern void X(codelet_hc2cfdftv_10)(planner *);
extern void X(codelet_hc2cfdftv_12)(planner *);
extern void X(codelet_hc2cfdftv_16)(planner *);
extern void X(codelet_hc2cfdftv_32)(planner *);
extern void X(codelet_hc2cfdftv_20)(planner *);
extern void X(codelet_hc2cbdftv_2)(planner *);
extern void X(codelet_hc2cbdftv_4)(planner *);
extern void X(codelet_hc2cbdftv_6)(planner *);
extern void X(codelet_hc2cbdftv_8)(planner *);
extern void X(codelet_hc2cbdftv_10)(planner *);
extern void X(codelet_hc2cbdftv_12)(planner *);
extern void X(codelet_hc2cbdftv_16)(planner *);
extern void X(codelet_hc2cbdftv_32)(planner *);
extern void X(codelet_hc2cbdftv_20)(planner *);


extern const solvtab X(solvtab_rdft_simd);
const solvtab X(solvtab_rdft_simd) = {
   SOLVTAB(X(codelet_hc2cfdftv_2)),
   SOLVTAB(X(codelet_hc2cfdftv_4)),
   SOLVTAB(X(codelet_hc2cfdftv_6)),
   SOLVTAB(X(codelet_hc2cfdftv_8)),
   SOLVTAB(X(codelet_hc2cfdftv_10)),
   SOLVTAB(X(codelet_hc2cfdftv_12)),
   SOLVTAB(X(codelet_hc2cfdftv_16)),
   SOLVTAB(X(codelet_hc2cfdftv_32)),
   SOLVTAB(X(codelet_hc2cfdftv_20)),
   SOLVTAB(X(codelet_hc2cbdftv_2)),
   SOLVTAB(X(codelet_hc2cbdftv_4)),
   SOLVTAB(X(codelet_hc2cbdftv_6)),
   SOLVTAB(X(codelet_hc2cbdftv_8)),
   SOLVTAB(X(codelet_hc2cbdftv_10)),
   SOLVTAB(X(codelet_hc2cbdftv_12)),
   SOLVTAB(X(codelet_hc2cbdftv_16)),
   SOLVTAB(X(codelet_hc2cbdftv_32)),
   SOLVTAB(X(codelet_hc2cbdftv_20)),
   SOLVTAB_END
};
