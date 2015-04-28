#!/usr/bin/env ruby
# In the following examples, the guide is: CGTTCTTATGGAAGCCGGGA followed by AGG PAM.
#
# In all these cases, take the genome position and add the matched length from CIGAR string to get the genome position of the mutation.
#
# Insertions only appear in the cigar string; but the inserted base can only be found from the read  e.g.
#DKNQZ:00014:00145       0       5       112767204       37      58M1I39M        *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCCGGGAAGGATCTGTATCAAGCCGTTCTGGAGAGTGCAG      CBDC??<;D1;?=?C@CC;;;;/;;;;;@C;;;;;-;;444*4@79@.22499:<@<:::/=D<D@D9=;;;;;;;44.2274;)-.-)--332226=      XT:A:U  NM:i:1  X0:i:1  X1:i:0  XM:i:0  XO:i:1  XG:i:1  MD:Z:97
# in the above case the inserted bases would be sequence_string[total matched length .. total matched length+insertion length] i.e. sequence_string[58..59] => "C"
#
# Deletions will be in both - they are ^N in the MD:Z tag, e.g
# DKNQZ:00025:00303       0       5       112767204       37      60M1D7M2I6M     *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA     CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112     XT:A:U  NM:i:3  X0:i:1  X1:i:0  XM:i:3  XO:i:1  XG:i:1  MD:Z:60^G13
# To get the base, need to parse the MD:Z tag and take the [ATGC] characters after the ^.
#
# In the case of a mismatch, the CIGAR string will be M (As this is made for alignment, not variant calling), whereas the reference base will appear in the MD:Z tag - e.g. here when CCgGG in the reference becomes CCCGG:
# DKNQZ:00167:00152       0       5       112767204       25      157M    *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCCGGAAGGATCTGTATCAAGCCGTTCTGGAGAGTGCAGTCCTGTTCCTATGGGTTCATTTCCAAGAAGAGGGTTTGTAAATGGAAGCAGAGGCTGAGG   ??CCCCACD>CC>CC>BB;///(00/--3--2222*222::=:C65882>>A::5:5:AC;;6:66).-33322529.29/22/9:2A522CC;:?CCEACCD@@@@DCDG@FA??<<;?@DDFDCC@@@=@D<@@CC>CC>?>>AC@;>6;A@;;7   XT:A:U  NM:i:7  X0:i:1  X1:i:0  XM:i:7  XO:i:0  XG:i:0  MD:Z:60G89A0A0A1T0A0C0
#
# ===>  Need the MD, CIGAR and reference sequence
#
# Methods for MD parsing:
# 	-cumulative reference length (eq. to matched_length in CIGAR).
# 	-organise into deletions & mismatches [no insertions in MD:Z]
# 	
require 'test/unit'
require '/Users/spettitt/work/scripts/crisprion/sam'
# DKNQZ:00025:00303       0       5       112767204       37      60M1D7M2I6M     *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA     CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112     XT:A:U  NM:i:3  X0:i:1  X1:i:0  XM:i:3  XO:i:1  XG:i:1  MD:Z:60^G13
