#!/bin/bash
# THIS BATCH RUNS MSA-RUN PROGRAM FROM IMAGIC version 2016-04-18

# ---------------- Parameters ----------------

msa_distance="MODULATION"                           # type of inter-image distances used
num_factors="25"                                    # number of eigenimages
num_iter="25"                                       # number of iterations
overcorrectionFactor="0.8"                          # overcorrection factor [0.0-0.9]
mpi_procs="1"                                       # if > 1, we run mpi version

# ------------------ Inputs ------------------

particles="input_particles"                         # input stack of aligned particles (should be at least centered)
mask="tmp/mask"                                         # mask file name

# ------------------ Outputs ------------------

msa_dir="MSA"                                       # MSA results are in this folder
eigen_img="eigen_img"                               # eigenimages
msa_pixvec_coord="msa_pixvec_coord"                 # eigenvectors in pixel space
msa_eigen_pixel="msa_eigen_pixel"                   # eigenvectors in image space

# -------------- END BATCH HEADER --------------

# check if we run mpi version
if [ ${mpi_procs} -eq "1" ]; then
    ${IMAGIC_ROOT}/msa/msa.e <<EOF
FRESH_MSA
${msa_distance}
${particles}
${mask}
${eigen_img}
NO
${num_iter}
${num_factors}
1
${overcorrectionFactor}
${msa_dir}/msa
NO
EOF

else
    ${IMAGIC_ROOT}/openmpi/bin/mpirun -np ${mpi_procs} -x IMAGIC_BATCH ${IMAGIC_ROOT}/msa/msa.e_mpi <<EOF
FRESH_MSA
${msa_distance}
${particles}
${mask}
${eigen_img}
NO
${num_iter}
${num_factors}
1
${overcorrectionFactor}
${msa_dir}/msa
YES
${mpi_procs}
NO_LOCAL_FILES
EOF

fi

# Modified 2016-04-25
#    2016-04-25 (gs) -- updated for new IMAGIC release, add mpi support
#    2016-02-26 (gs) -- first version
