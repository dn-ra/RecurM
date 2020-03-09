#!/bin/bash
#edits contig names to include the assembly that they are sourced from
#for use in nucmer for repeatm

#remember that the assembly name is delimited from the node by __ (double underscore) this is how RecurM associates contigs with an assembly

#input: list of assemblies 
#output: separate fasta file for each assembly (filtered to your desired length), names appended, with contig names containing node-assembly links (eg. >ASSEMBLY__CONTIG))

module load parallel

function amend_names {
assembly=$(basename $1)
echo $1 $assembly
sed  's/>/&'"$assembly"'__/' $1 > node_assembly_linked."$assembly" #stored in current directory
echo $1 amended
}
export -f amend_names


ls $@ > source_assemblies.txt

parallel --env amend_names -a source_assemblies.txt --verbose amend_names {}

wait


module unload parallel
