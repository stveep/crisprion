require 'bio'
require '/Users/spettitt/work/scripts/cigar/iterate_pairs'
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
		h = Hash.new{|k,v| h[k] = ""}
		h[:unsplit_tags] = @tags.join(" ")
		# Reconstruct the tag and assign a value. Assumes all tags are of the form X:Y:stuff
		@tags.map!{|tag| tag.split(/:/,3)}
		@tags.each{|a,b,c| h[a + ":" + b] = c}
		@tags = h
		@pos = @pos.to_i
		@mapq = @mapq.to_i
		@mrnm == "*" ? @mrnm = nil : @mrnm = @mrnm
		@mpos == "0" ? @mpos = nil : @mpos = @mpos.to_i
		@isize == "0" ? @isize = nil : @isize = @isize.to_i
		@seq = Bio::Sequence::NA.new(@seq)
	end

end

class Bio::Alignment::SAM::MDZ
	include Bio::Alignment::IteratePairs
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

	# Sums the total length of the reference sequence represented by the MD:Z tag
	def ref_length
		count_length(@tag)	
	end

	def subregion(offset,length)
		@breakdown = @tag.scan(/\^[ATGC]+|\d+|[GATC]/) #=> array of number, deletion or base elements
		@breakdown.map!{|a| [a,count_length(a)]}
		# Breakdown is now an array of [operation, length] arrays
		# Use IteratePairs mixin (from CIGAR parser)
		new_pairs = iterate_pairs(@breakdown,offset,length)
		# Need to also modify operator part of pair in the case of matches
		new_pairs.map!{|a| if a[0].match(/^\d+$/) then [a[1],a[1]] else a end }	
		new_pairs.map!{|a| a[0]}.join("")
	end

	
	
	private 
	def count_length(tag)
		#Need the sum of all "movement" operations (i.e. numbers) as well as any substituted bases (count 1 each)
		if tag =~ /^\d+$/
			tag.to_i
		else
			temp_tag = tag.dup
			temp_tag.gsub!(/\^/,"")  # Deletions need to be counted - sub the caret character out and count the remaining base characters
			movements = temp_tag.split(/[GATC]+/).map(&:to_i).reduce(:+) # Sum numbers
			deletions = temp_tag.split(/\d+/).map(&:length).reduce(:+) # Sum number of base chars
			# If either one of these is nil, to_i coerces to zero
			movements.to_i + deletions.to_i 
		end
	end

end


