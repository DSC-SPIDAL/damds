#!/usr/bin/env bash

x='x'
#opts="-XX:+UseConcMarkSweepGC -XX:ParallelCMSThreads=4 -Xms2G -Xmx2G"
opts="-XX:+UseG1GC -Xms512m -Xmx512m"
tpn=1
wd=`pwd`
summary=summary.txt
timing=timing.txt
p=$9
$BUILD/bin/mpirun --report-bindings --mca btl ^tcp -n $p --hostfile nodes.txt java $opts -cp ../target/damds-1.1-jar-with-dependencies.jar -DDistanceMatrixFile=$1 -DWeightMatrixFile=$2 -DPointsFile=$3.txt -DNumberDataPoints=$4 -DTimingFile=$5 -DSummaryFile=$6 -DInitialPointsFile=$7 edu.indiana.soic.spidal.damds.Program -c $8 -n $p -t $tpn | tee $4.summary.txt
echo "Finished $0 on `date`" >> status.txt

