#!/bin/bash

for SPECIES in Hanseniaspora_uvarum	Lachancea_fermentati Hanseniaspora_vineae Lachancea_kluyveri Kluyveromyces_lactis Lachancea_thermotolerans Kluyveromyces_marxianus
do


mafft --auto $SPECIES\.fasta > ../alnmnts/$SPECIES\-aln.fasta

### raxML gene tree
~/standard-RAxML/raxmlHPC-SSE3 -f a -m GTRCAT -p 12345 -x 12345 -# 100 -s ../alnmnts/$SPECIES\-aln.fasta -n $SPECIES

### Using Newick utilities to extract branch lengths
mv RAxML_bestTree* ../trees/$SPECIES\.tree

rm RAxML*

done

