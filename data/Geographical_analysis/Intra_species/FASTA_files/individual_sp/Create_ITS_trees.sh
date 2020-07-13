#!/bin/bash

for SPECIES in Aureobasidium_pullulans Candida_boidinii Candida_pseudolambica Candida_railenensis Candida_sake Candida_santamariae Candida_sp Candida_vartiovaarae Create_ITS_trees.sh Cryptococcus_flavescens Cryptococcus_gastricus Cryptococcus_pinus Cryptococcus_sp Curvibasidium_cygneicollum Cyberlindnera_saturnus Cyberlindnera_veronae Cystofilobasidium_capitatum Debaryomyces_hansenii Geotrichum_candidum Hanseniaspora_uvarum Hanseniaspora_vineae Kazachstania_martiniae Kluyveromyces_lactis Kluyveromyces_marxianus Kodamaea_meredithiae Lachancea_fermentati Lachancea_kluyveri Lachancea_thermotolerans Lecythophora_sp Leucosporidium_scottii Metschnikowia_pulcherrima_sp._complex Metschnikowia_reukaufii Metschnikowia_viticola Meyerozyma_caribbica Mrakia_blollopis Mrakia_gelida Mrakia_sp Nakazawaea_anatomiae:populi_sp._complex Nakazawaea_holstii Papiliotrema_flavescens Pichia_fermentans Pichia_kudriavzevii Pichia_manshurica Rhodosporidiobolus_colostri Rhodotorula_nothofagi Saccharomyces_cerevisiae Saccharomyces_eubayanus Saccharomyces_paradoxus Scheffersomyces_ergatensis Schwanniomyces_pseudopolymorphus Teunomyces_carpophila Teunomyces_cretensis:kruisii_complex Torulaspora_delbrueckii Trichosporon_porosum Wickerhamomyces_anomalus
do

#below line removes sequences less than minimum of 200 bp
awk -v min="200" 'BEGIN {RS = ">" ; ORS = ""} length($2) >= min {print ">"$0}' $SPECIES\.fasta > $SPECIES\-200bpmin.fasta

#below lines remove sequences with more than 3 Ns
~/anaconda3/bin/faidx --transform nucleotide $SPECIES\-200bpmin.fasta > base_counts.txt

cat base_counts.txt | awk '{if ($8 < 3); print $1}' > seqs_tokeep.txt

grep -w -A 2 -f  seqs_tokeep.txt $SPECIES\-200bpmin > $SPECIES\-Nreduced.txt

sed '/^--$/d' $SPECIES\-Nreduced.txt > $SPECIES\-Nreduced.fasta

rm *Nreduced.txt
rm base_counts.txt
rm seqs_tokeep.txt

mafft --auto $SPECIES\-200bpmin.fasta > ../alnmnts/$SPECIES\-aln.fasta

#using TrimAl gappyout method to trim sequences
~/anaconda3/bin/trimal -in ../alnmnts/$SPECIES\-aln.fasta -out ../alnmnts/$SPECIES\-aln-reduced.fasta -gappyout


### raxML gene tree
~/standard-RAxML/raxmlHPC-SSE3 -f a -m GTRCAT -p 12345 -x 12345 -# 100 -s ../alnmnts/$SPECIES\-aln-reduced.fasta -n $SPECIES

### Using Newick utilities to extract branch lengths
mv RAxML_bestTree* ../trees/$SPECIES\-ML.tree
mv RAxML_bipartitions.* ../trees/$SPECIES\-ML-BS.tree

rm RAxML*

done

