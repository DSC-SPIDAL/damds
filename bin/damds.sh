#!/usr/bin/env bash

x='x'
#opts="-XX:+UseConcMarkSweepGC -XX:ParallelCMSThreads=4 -Xms2G -Xmx2G"
opts="-XX:+UseG1GC -Xms512m -Xmx512m"
tpn=1
wd=`pwd`

$BUILD/bin/mpirun --report-bindings --mca btl ^tcp -n 1 --hostfile nodes.txt java $opts -cp ../target/damds-1.1-jar-with-dependencies.jar -DNumberDataPoints=$2 -DDistanceMatrixFile=$1 -DPointsFile=$3.txt -DWeightMatrixFile=$4 -DTimingFile=timing.txt -DSummaryFile=summary.txt -DInitialPointsFile=init_points_37.txt edu.indiana.soic.spidal.damds.Program -c config.properties -n 1 -t 1 | tee $4.summary.txt
echo "Finished $0 on `date`" >> status.txt

