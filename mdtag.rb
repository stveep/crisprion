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
require 'rubygems'
require 'bio'

class Bio::Alignment::SAM
	# Field names from SAM specification
	attr_accessor :qname, :flag, :rname, :pos, :mapq, :cigar, :mrnm, :mpos, :isize, :seq, :qual, :tags
	# Aliases included for more intuitive naming when using for genomic alignments:
	alias_method :chr, :rname
	alias_method :opt, :tags # vice versa as opt is the "proper" sam name for the tag fields

	def initialize(line)
		@fields = line.split(/\s+/)
		@qname, @flag, @rname, @pos, @mapq, @cigar, @mrnm, @mpos, @isize, @seq, @qual = @fields
		@tags = @fields[11..-1]
		@h = Hash.new{|k,v| @h[k] = ""}
		@h[:unsplit_tags] = @tags.join(" ")
		# Reconstruct the tag and assign a value. Assumes all tags are of the form X:Y:stuff
		@tags.map!{|tag| tag.split(/:/,3)}
		@tags.each{|a,b,c| @h[a + ":" + b] = c}
		@tags = @h
		@pos = @pos.to_i
		@mapq = @mapq.to_i
		@mrnm == "*" ? @mrnm = nil : @mrnm = @mrnm
		@mpos == "0" ? @mpos = nil : @mpos = @mpos.to_i
		@isize == "0" ? @isize = nil : @isize = @isize.to_i
		@seq = Bio::Sequence::NA.new(@seq)
	end

end

class Bio::Alignment::SAM::MDZ
	attr_accessor :tag, :variants	
	@@regexp = /MD:Z:([\w^]+)/
	@@format = /[\w^]+/
	def initialize(string)
		if string.match(@@regexp)
			@tag = $~[1]
		elsif string.match(@@format)
			#Assume tag given without MD:Z: leader
			@tag = string
		else
			raise "Tag not of expected format."
		end
	end

	# Sums the total length of the reference sequence represented by the MD:Z tag (or part of)
	def ref_length
		#Need the sum of all "movement" operations (i.e. numbers) as well as any substituted bases (count 1 each)
		if @tag =~ /^\d+$/
			@tag.to_i
		else
			temp_tag = @tag.dup
			temp_tag.gsub!(/\^/,"")  # Deletions need to be counted - sub the caret character out and count the remaining base characters
			movements = temp_tag.split(/[GATC]+/).map(&:to_i).reduce(:+) # Sum numbers
			deletions = temp_tag.split(/\d+/).map(&:length).reduce(:+) # Sum number of base chars
			movements + deletions 
		end
	end

end


# DKNQZ:00025:00303       0       5       112767204       37      60M1D7M2I6M     *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA     CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112     XT:A:U  NM:i:3  X0:i:1  X1:i:0  XM:i:3  XO:i:1  XG:i:1  MD:Z:60^G13
class SAMTest < Test::Unit::TestCase
	def test_split_types
		sam = Bio::Alignment::SAM.new("DKNQZ:00025:00303       0       5       112767204       37      60M1D7M2I6M     *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA     CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112     XT:A:U NM:i:3 X0:i:1 X1:i:0 XM:i:3 XO:i:1 XG:i:1 MD:Z:60^G13")
		assert(sam.qname.is_a? String)
		assert(sam.flag.is_a? String)
		assert(sam.rname.is_a? String)
		assert(sam.pos.is_a? Integer)
		assert(sam.mapq.is_a? Integer)
		assert(sam.cigar.is_a? String)
		assert(sam.mrnm.is_a?(String) || sam.mrnm.nil?)
		assert(sam.mpos.is_a?(Integer) || sam.mpos.nil?)
		assert(sam.isize.is_a?(Integer) || sam.isize.nil?)
		assert(sam.seq.instance_of? Bio::Sequence::NA)
		assert(sam.qual.is_a? String)
		assert(sam.tags.is_a? Hash)

	end
	def test_split
		sam = Bio::Alignment::SAM.new("DKNQZ:00025:00303       0       5       112767204       37      60M1D7M2I6M     *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA     CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112     XT:A:U  NM:i:3  X0:i:1  X1:i:0  XM:i:3  XO:i:1  XG:i:1  MD:Z:60^G13")
		assert_equal(sam.qname,"DKNQZ:00025:00303","ID not as expected")
		assert_equal(sam.flag,"0","Flag not as expected")
		assert_equal(sam.rname,"5","Chr not as expected")
		assert_equal(sam.pos,112767204,"Position not as expected")
		assert_equal(sam.mapq,37,"Quality not as expected")
		assert_equal(sam.cigar,"60M1D7M2I6M","CIGAR not as expected")
		assert_equal(sam.mrnm,nil,"Mate chr not as expected [* -> nil]")
		assert_equal(sam.mpos,nil,"Mate pos not as expected [0 -> nil]")
		assert_equal(sam.isize,nil,"Insertion size not as expected [0 -> nil]")
		assert_equal(sam.seq,Bio::Sequence::NA.new("GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA"),"Sequence not as expected")
		assert_equal(sam.qual,"CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112","Base quality string not as expected")
 assert_equal(sam.tags,{:unsplit_tags => "XT:A:U NM:i:3 X0:i:1 X1:i:0 XM:i:3 XO:i:1 XG:i:1 MD:Z:60^G13", "XT:A" => "U", "NM:i" => "3", "X0:i" => "1", "X1:i" => "0", "XM:i" => "3", "XO:i" => "1","XG:i" => "1", "MD:Z" =>"60^G13"})
		
	end
	
	def test_aliases
		sam = Bio::Alignment::SAM.new("DKNQZ:00025:00303       0       5       112767204       37      60M1D7M2I6M     *       0       0       GCAGTAATTTCCCTGGAGTAAAACTGCGGTCAAAAATGTCCCTCCGTTCTTATGGAAGCCGGAAGGAAGTCTGTA     CCCCCC@CE>CC<CC@CB;;;;.;;;;;AC;::::+:92A:=CCAEE=?>;=:@<B?:<6<*/*/*/*/911112     XT:A:U  NM:i:3  X0:i:1  X1:i:0  XM:i:3  XO:i:1  XG:i:1  MD:Z:60^G13")
		assert_equal(sam.chr,sam.rname)
		assert_equal(sam.tags,sam.opt)
	end

end


class MDZTest < Test::Unit::TestCase
	def test_match
		mdz = Bio::Alignment::SAM::MDZ.new("MD:Z:60^G13")
		assert_equal mdz.tag, "60^G13"
	end

	def test_ref_length
		mdz_del = Bio::Alignment::SAM::MDZ.new("MD:Z:60^G13")
		assert_equal(74, mdz_del.ref_length)	
		mdz_sub = Bio::Alignment::SAM::MDZ.new("MD:Z:60G0A13")	
		assert_equal(75, mdz_sub.ref_length)	
		mdz_none = Bio::Alignment::SAM::MDZ.new("MD:Z:150")
		assert_equal(150, mdz_none.ref_length)	

	end

end
