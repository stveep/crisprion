require "/home/breakthr/spettitt/scratch/scripts/fastq.rb"
require "open4"
fastq_files = Dir.glob(ARGV)
working_directory = Dir.getwd
fastq_files = fastq_files.map!{|a| working_directory + "/" + File.basename(a) }
fastq_stems = fastq_files.map{|a| a.sub(/\.fastq$/,"") }

use_lsf = true
project = ENV["LSF_PROJECT"]
start_clip = /\w+AGACTATCTTTC/
end_clip = /not used/

min_length = 50
bwa = "/apps/bwa/0.7.5a/bwa"
samtools = "/apps/samtools/1.1/bin/samtools"
#bwaindex ="/home/breakthr/spettitt/scratch/genome/human/chromosomes/Homo_sapiens.GRCh38.dna_rm.chromosome.5.fa"
bwaindex ="/home/breakthr/spettitt/scratch/genome/hsvTK/TNP.fa"
region = "NULL"

desc "Remove primer sequence from start of fastq reads"
task :clip => fastq_stems.map{|i| i+".fastqclip"}

desc "Produce .sai files (1st step of alignment)"
task :align1 => fastq_stems.map{|i| i+".sai"}

desc "Produce .sam files (2nd step of alignment)"
task :align2 => fastq_stems.map{|i| i+".bam"}

desc "Sort bam file"
task :sort => fastq_stems.map{|i| i+".sort"}

desc "Pileup"
task :bcf => fastq_stems.map{|i| i+".bcf"}

desc "Generate sam files from bam files"
task :sam => fastq_stems.map{|i| i+".sam"}

desc "Cleanup temporary LSF files (.otmp and .etmp)"
task :clean do |t|
	outfiles = Dir.glob("*.otmp") + Dir.glob("*.etmp")
	outfiles.each do |file|
		File.delete(file)	
	end
end

rule('.fastqclip' => '.fastq') do |t|
	f = Fastq.new(t.prerequisites[0])
	outfile = File.new(t.name,"w")	
	f.each do |r|
		if (r.seq =~ start_clip)
		clip = start_clip.match(r.seq).end(0)
		# Any adapter at the 3' end? 
		eclip = -1
		if r.seq.match(end_clip)
			#offset(0)[1] gives the end position of the first match
			eclip = $~.offset(0)[1]-1
		end
		next if (eclip > 0 && (eclip - clip < 1))
		next if r.seq[clip..eclip].length != r.qual[clip..eclip].length
		outfile.puts(r.name.chomp + "\n" + r.seq[clip..eclip] + "\n+\n" + r.qual[clip..eclip])
		end

	end

end

rule('.sai' => '.fastqclip') do |t|
        bwa_command = "'"+"#{bwa} aln #{bwaindex} #{t.prerequisites[0]} > #{t.name}"+"'"
	if use_lsf
		system("bsub -P #{project} -J #{t.name} -o #{t.name+".otmp"} -e #{t.name+".etmp"} #{bwa_command}")
	else
		system(bwa_command)
	end
end

rule('.bam' => [ proc {|taskname| taskname.sub('.bam','.sai')}\
                ,proc {|taskname| taskname.sub('.bam','.fastqclip')}]) do |t|
	bwa_command = "'"+"#{bwa} samse #{bwaindex} #{t.prerequisites[0]} #{t.prerequisites[1]} | #{samtools} view -o #{t.name} -bS -F 0x04 - "+"'"
	if use_lsf
		system("bsub -P #{project} -J #{t.name} -o #{t.name+".otmp"} -e #{t.name+".etmp"} #{bwa_command}")
	else
		system(bwa_command)
	end
end

rule('.sort' => '.bam') do |t|
	command = "#{samtools} sort -f #{t.prerequisites[0]} #{t.name}"
	if use_lsf
		system("bsub -P #{project} -J #{t.name} -o #{t.name+".otmp"} -e #{t.name+".etmp"} #{command}")
	else
		system(command)
	end
end

rule('.sam' => '.bam') do |t|
	command = "#{samtools} view -m #{min_length} #{t.prerequisites[0]} > #{t.name}"
	pid, stdin, stout, sterr = Open4::popen4 command 
	Process::waitpid2 pid
end 

rule('.bcf' => '.sort') do |t|
	index_command = "#{samtools} index #{t.prerequisites[0]}"	
	pid, stdin, stout, sterr = Open4::popen4 index_command 
	# Wait for the process to finish - need the index file to be there before starting pileup
	status = Process::waitpid2 pid
	command = "'"+"#{samtools} mpileup -gu -f #{bwaindex} -r #{region} #{t.prerequisites[0]} > #{t.name}"+"'"
	if use_lsf
		system("bsub -P #{project} -J #{t.name} -o #{t.name+".otmp"} -e #{t.name+".etmp"} #{command}")
	else
		system(command)
	end
end
