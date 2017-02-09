#!/usr/bin/env bash

jar=damds-1.1-jar-with-dependencies.jar
x='x'
opts="-XX:+UseG1GC -Xms768m -Xmx1024m"
tpn=1
wd=`pwd`
summary=summary.txt
timing=timing.txt
$BUILD/bin/mpirun --report-bindings --mca btl ^tcp -np $3 --hostfile nodes.txt java $opts -cp ../target/$jar edu.indiana.soic.spidal.damds.Program -c $1 -n $2 -t 1 | tee summary.txt
echo "Finished $0 on `date`" >> status.txt
