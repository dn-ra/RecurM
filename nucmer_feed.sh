#!/bin/bash

#This script performs nucmer alignment of all assemblies in input
#input assemblies must be node-assembly matched (created by "nucmer_ready_assemblies.sh")

#input: list of n assemblies to be processed
#ouput: n-1 .delta files to be processed by RecurM

echo Script initiated at $(date)

module load mummer/4.0.0beta2

read -a ASSEMBLIES <<< ${@}

echo $# assemblies submited for alignment

mkdir nucmer_feed_out
touch assemblies_iter_built.fa #this is the query sequences file. Increases by adding the previous reference file for each iteration.

QUERY_SIZE=0

for ((i=0; i<$#-1; i++)); do
CURR_QUE="${ASSEMBLIES["$i"]}";
CURR_REF="${ASSEMBLIES[$(expr $i + 1)]}";
ITER=$(expr $i + 1);

SEQ_COUNT=$(grep '>' $CURR_QUE | wc -l);
QUERY_SIZE=$(expr $QUERY_SIZE + $SEQ_COUNT);

echo adding "$SEQ_COUNT" to concatenated query file in iteration "$ITER";
cat $CURR_QUE >> assemblies_iter_built.fa;


echo processing "$QUERY_SIZE" sequences in query file against reference file "$CURR_REF";

/usr/bin/time nucmer "$CURR_REF" assemblies_iter_built.fa --delta=nucmer_feed_out/"$(date +%d%m%y)"_nucmer_iter_"$ITER".delta -t 24;
done

rm assemblies_iter_built.fa
