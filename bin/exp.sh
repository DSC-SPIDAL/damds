#!/usr/bin/env bash

INPUT_DIR=$1
OUT_DIR=$2
CONF_FILE=/N/u/skamburu/data/flink/mpi_full/config.properties
P=$3
N=$4
#points=( 1000 2000 4000 8000 16000 32000 )
points=( 1000 )

for i in "${points[@]}"
do
    mkdir -p $OUT_DIR/$i
    #sh damds.sh $INPUT_DIR/$i/dist.bin $INPUT_DIR/$i/weight.bin $OUT_DIR/$i/points $i $OUT_DIR/$i/timing.txt $OUT_DIR/$i/summary.txt $INPUT_DIR/$i/init_points $CONF_FILE $P $N
    sh damds.sh /N/u/skamburu/projects/damds_lrt/bin/whiten_dataonly_fullData.2320738e24c8993d6d723d2fd05ca2393bb25fa4.4mer.dist.c#_1000.bin "" $OUT_DIR/$i/points $i $OUT_DIR/$i/timing.txt $OUT_DIR/$i/summary.txt "" $CONF_FILE $P $N
done
