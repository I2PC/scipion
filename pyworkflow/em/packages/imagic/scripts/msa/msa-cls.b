#!/bin/bash
# THIS BATCH RUNS MSA-CL AND MSA-SUM PROGRAMS FROM IMAGIC version 2016-04-18

# ---------------- Parameters ----------------

eigs_num="15"                           # number of eigenimages used
cls_num="10"                            # number of classes
perc_ign="0"                            # percent of images to ignore
downweight="NO"                         # downweight small classes?
perc_ign_bad="0"                        # percent of worst class members to ignore
classno="108"                           # header index CLASSNO (Class number in MSA classification )

# ------------------ Inputs ------------------

particles="tmp/particles"             # input stack particles after MSA-RUN

# ------------------ Outputs ------------------

cls_dir="MSA-cls"                       # results are in this folder
msa_cls_img="classes"                   # output rootname


# -------------- END BATCH HEADER --------------

echo "IMAGIC program: msa-cl (v. 11-Feb-2016) is starting"
${IMAGIC_ROOT}/msa/classify.e <<EOF
IMAGES/VOLUMES
HAC
${particles}
${perc_ign}
${eigs_num}
YES
${cls_num}
${cls_dir}/${msa_cls_img}
EOF

echo "IMAGIC program: msa-sum (v. 19-Aug-2015) is starting"
${IMAGIC_ROOT}/msa/classum.e <<EOF
${particles}
${cls_dir}/${msa_cls_img}
${cls_dir}/${msa_cls_img}_avg
${downweight}
${perc_ign_bad}
NONE
EOF

echo "IMAGIC program: headers (v. 21-Jan-2016) is starting"
${IMAGIC_ROOT}/stand/headers.e <<EOF > /dev/null
PLT_OUT
INDEX/LABEL
NUMBER_OF_INDEX
${classno}
YES
${particles}
${cls_dir}/class_assignment.plt
EOF


# Modified 2016-04-25
#    2016-04-25 (gs) -- updated for new IMAGIC release
#    2016-03-04 (gs) -- added printing of MSA class numbers
#    2016-03-03 (gs) -- first version
