#!/bin/bash

source setup.sh
#echo ========== TEST ==========
#echo LD_LIBRARY_PATH
#echo $LD_LIBRARY_PATH
#echo "ls /usr/local/jetreader_build/lib/libyaml*"
#ls /usr/local/jetreader_build/lib/libyaml*
#echo INPUTFILECOUNT
#echo $INPUTFILECOUNT
#echo FILELIST
#echo $FILELIST
#cp $FILELIST /gpfs01/star/pwg/youqi/run15/.
#echo ========== TEST ==========
./bin/jets -1 $INPUTFILECOUNT $FILELIST $JOBID
#./bin/jets -1 $INPUTFILECOUNT "/tmp/youqi/$JOBID/INPUTFILES/st_physics_16064077_raw_5500002.picoDst.root" "/gpfs01/star/pwg/youqi/run15/result/" $JOBID # this works
#./bin/jets -1 $INPUTFILECOUNT "st_physics_16064077_raw_5500002.picoDst.root" "/gpfs01/star/pwg/youqi/run15/result/" $JOBID # this works 
