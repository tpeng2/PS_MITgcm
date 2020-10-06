#!/bin/sh
GPNAME=$1
CNAME=$2

dtpath=/home/tpeng/scratch/MITgcm/examples/${GPNAME}/${CNAME}/
oppath=/home/tpeng/MITgcm/postproc/results/${GPNAME}/${CNAME}/
bkpath=/home/tpeng/MITgcm/postproc/results/${GPNAME}/${CNAME}/backup/

# create folder if it doesn't exist
mkdir -p ${oppath}
mkdir -p ${bkpath}
# copy files
echo 'Copy complete grids'
sleep 2
cp ${dtpath}XC.*ta ${oppath}
cp ${dtpath}YC.*ta ${oppath}
cp ${dtpath}XG.*ta ${oppath}
cp ${dtpath}YG.*ta ${oppath}
cp ${dtpath}RC.*ta ${oppath}
cp ${dtpath}RF.*ta ${oppath}

echo 'Copy incremental grids'
sleep 2

cp ${dtpath}D*C.*ta ${oppath}
cp ${dtpath}D*G.*ta ${oppath}
cp ${dtpath}D*F.*ta ${oppath}

echo 'Copy hFac files'
sleep 2
cp ${dtpath}hFac*ta ${oppath}

echo 'Copy forcing and topography'
sleep 2
cp ${dtpath}*.box ${oppath}
cp ${dtpath}*.sin_y ${oppath}

# back up data
echo 'Back up configuration'
sleep 2
cp ${dtpath}data* ${bkpath}
cp ${dtpath}eedata ${bkpath}
cp ${dtpath}*.box ${bkpath}
