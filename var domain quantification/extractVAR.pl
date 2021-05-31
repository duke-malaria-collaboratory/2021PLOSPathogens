#!/usr/bin/perl
##
##
##use warnings;

use strict;

if ($#ARGV != 4) {
    print "usage: Enter full path to target genome directory, output directory, directory with R1 reads, directory with R2 reads, and GTF file \ne.g. countRNASeqReads.pl /home/genome /home/out /home/reads/1 /home/reads/2 /home/file.gtf";
    exit;
}
my $GenomeDir = $ARGV[0];
my $out = $ARGV[1];
my $pair1 = $ARGV[2];
my $pair2 = $ARGV[3];
my $GTF = $ARGV[4];

chomp $GenomeDir;
if ( -d $GenomeDir )
{
    print "\nGenome Directory: $GenomeDir\n";
}
else { die "\nGenome Directory does not exist\n"; }

getReference("$GenomeDir");

chomp $out;
makeOutDirs("$out");

chomp $pair1;
chomp $pair2;
chomp $GTF;

my @PairedReadfiles1 = ();
my @PairedReadfiles2 = ();
my $trim1;
my $trim2;
my @trimReads1 = ();
my @trimReads2 = ();


@PairedReadfiles1 = getReads("$pair1");
@PairedReadfiles2 = getReads("$pair2");

align($GenomeDir, "$pair1", "$pair2", "$out", "$GTF", \@PairedReadfiles1, \@PairedReadfiles2);

############################################################
sub getReference
{
    my $GenomeDir = shift();
    my @Files = ();
    my $file;
    my @Ext = ();
    my $fasta;
    my $ext;
    my $i = 0;
    my $found = 0;
    my $GTFfile;
    my $Index;
    if (-d $GenomeDir)
    {
        opendir (GENOMEDIR, "$GenomeDir");
        @Files = readdir GENOMEDIR;
        closedir GENOMEDIR;
        splice (@Files, 0, 2);
        foreach $file (@Files)
        {
            ($Ext[$i]) = $file =~ /(\.[^.]+)$/;
			if ($Ext[$i] eq "")
			{
				$i++;
			}
            elsif ((($Ext[$i] eq ".fasta") || ($Ext[$i] eq ".fa")) && ($found == 0))
            	{
                	$fasta = $file;
                	$found = 1;
                	print "\n\tUsing genome FASTA file $GenomeDir/$fasta\n";
            	}
            else{$i++};
        }
    }
    else {die "\n\t**Reference directory not found**\n";}

}
############################################################
sub makeOutDirs
{
    my $topOut = shift();
    unless (-d $topOut)
    {system "mkdir $topOut";}
    if ( -d $topOut)
    {
        unless (-d "$topOut/fastqc_in")
        {system "mkdir $topOut/fastqc_in";}
        print "\n\tFastqc .html files of input reads stored in $topOut/fastqc_in\n";

        unless (-d "$topOut/STAR_out")
        {system "mkdir $topOut/STAR_out";}
        print "\tTemporary .bam files stored in $topOut/STAR_out\n";

        unless (-d "$topOut/trim")
        {system "mkdir $topOut/trim";
        system "mkdir $topOut/trim/1";
        system "mkdir $topOut/trim/2";
        system "mkdir $topOut/trim/singleton";
        system "mkdir $topOut/trim/singleton/1";
        system "mkdir $topOut/trim/singleton/2";
        system "mkdir $topOut/trim/Log";
        system "mkdir $topOut/trim/Summary";
        system "mkdir $topOut/trim/fastqcTrim";}
        print "\tTrimmed reads stored in $topOut/trim\n";

        unless (-d "$topOut/fastqc_trim")
        {system "mkdir $topOut/fastqc_trim";}
        print "\n\tFastqc .html files of trimmed reads stored in $topOut/fastqc_trim\n";
        
        unless (-d "$topOut/quant")
        {system "mkdir $topOut/quant";}
        print "\tSalmon quant files stored in $topOut/quant\n";
        
        unless (-d "$topOut/sort_bam")
        {system "mkdir $topOut/sort_bam";}
        print "\tSorted and indexed .bam files stored in $topOut/sort_bam\n";
        
        unless (-d "$topOut/varReads")
        {system "mkdir $topOut/varReads";
        system "mkdir $topOut/varReads/bam";
        system "mkdir $topOut/varReads/1";
        system "mkdir $topOut/varReads/2";
        system "mkdir $topOut/varReads/singleton";}
        print "\tVAR reads stored in $topOut/varReads\n";

        unless (-d "$topOut/UnmappedReads")
        {system "mkdir $topOut/UnmappedReads";
	system "mkdir $topOut/UnmappedReads/1";
	system "mkdir $topOut/UnmappedReads/2";}
        print "\tUnmapped reads stored in $topOut/UnmappedReads\n";

    }
}
############################################################
sub getReads
{
    my $readsDir = shift();
    my @Reads = ();
    my $elem;
    if ( -d $readsDir)
    {
        opendir (READS, "$readsDir");
        @Reads = readdir READS;
        closedir READS;
        splice (@Reads, 0, 2);
    }
    else { die "\n\t**Reads not found**\n"; }

    my @sortReads = sort(@Reads);
    return @sortReads;
    }
############################################################
sub QCreads
{
    my $read1Dir = shift;
    my $read2Dir = shift;
    my $outDir = shift;
    my $pair1Reads = shift;
    my $pair2Reads = shift;
    my $format = shift;
    my @Reads1 = @{$pair1Reads};
    my @Reads2 = @{$pair2Reads};
    my $i;
    my $size = @Reads1;
    my $size2 = @Reads2;
    if ($size != $size2) {die "Paired end read files unequal";}

    for ($i = 0; $i < $size; $i++)
    {
        system "fastqc -o $outDir/fastqc_in -f $format $read1Dir/$Reads1[$i] $read2Dir/$Reads2[$i]";
    }
}
############################################################
sub align
{
    my $GenomeDir = shift;
    my $read1Dir = shift;
    my $read2Dir = shift;
    my $out = shift;
    my $GTF = shift;
    my $pair1Reads = shift;
    my $pair2Reads = shift;
    my @Reads1 = @{$pair1Reads};
    my @Reads2 = @{$pair2Reads};
    my $i;
    my $j;
    my $size = @Reads1;
    my $size2 = @Reads2;
    my @PRE;
    my $prefix;
    my @bed = ("1:29510-37126", "1:42367-46507", "1:607390-614893", "2:25232-31168", "2:916352-923648", "3:36965-44482", "3:1030822-1038254", "4:28706-37677", "4:45555-56860", "4:545987-553810", "4:561667-569342", "4:576810-584668", "4:591949-599849", "4:935031-941875", "4:946169-953773", "4:958067-965611", "4:969031-976591", "4:1156423-1167821", "4:1173053-1180226", "5:20929-28456", "5:1333470-1342964", "6:3503-12835", "6:18586-22721", "6:723117-731490", "6:1353946-1366430", "6:1374797-1382627", "7:20307-28083", "7:511950-519425", "7:527338-535131", "7:542906-550520", "7:552158-559078", "7:566726-574308", "7:581386-588923", "7:590326-597733", "7:1417588-1426234", "8:21361-28653", "8:29700-39033", "8:40948-50939", "8:431165-439051", "8:440408-448062", "8:459312-466567", "8:1435794-1443449", "9:20080-27885", "9:1486070-1490187", "9:1495579-1503336", "10:28490-36164", "10:1642401-1649948", "11:24160-31598", "11:32666-42386", "11:2025817-2035886", "12:16973-24497", "12:32703-41940", "12:46788-56805", "12:766654-774197", "12:1694150-1703087", "12:1704512-1712490", "12:1719574-1727456", "12:1735543-1743408", "12:2241271-2248962", "13:21364-28787", "13:33959-44742", "13:60097-61172", "13:2884786-2892340");
	my $numVars = @bed;

    if ($size != $size2) {die "Paired-end read file names unequal";}

    for ($i = 0; $i < $size; $i++)
    {
        @PRE = split(/\./, $Reads1[$i]);
        $prefix = $PRE[0];

        print "\n$prefix\n";
        system "mkdir $out/quant/$prefix";

        my $read1 = "$read1Dir/$prefix" . ".1.fastq.gz";
        my $read2 = "$read2Dir/$prefix" . ".2.fastq.gz";
		
       	system "STAR --runMode alignReads --genomeDir $GenomeDir --readFilesIn $read1 $read2 --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --runThreadN 12 --alignIntronMin 5 --alignIntronMax 3000 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0  --outFilterMatchNmin 15 --outFilterMismatchNmax 10 --outReadsUnmapped Fastx --outFileNamePrefix $out/STAR_out/$prefix.";

	system "samtools sort -o $out/sort_bam/$prefix.sort.bam $out/STAR_out/$prefix.Aligned.out.bam";
		
	system "rm $out/STAR_out/$prefix.Aligned.out.bam";

	system "mv $out/STAR_out/$prefix.Unmapped.out.mate1 $out/UnmappedReads/1/$prefix.1.fastq";

	system "mv $out/STAR_out/$prefix.Unmapped.out.mate2 $out/UnmappedReads/2/$prefix.2.fastq";
		
	system "samtools index -b $out/sort_bam/$prefix.sort.bam";
	
	system "touch $out/varReads/bam/bams_$prefix.txt";	
	for ($j = 0; $j < $numVars; $j++)
		{system "samtools view -b $out/sort_bam/$prefix.sort.bam $bed[$j] > $out/varReads/bam/var_$prefix.$j.bam";
		system "echo $out/varReads/bam/var_$prefix.$j.bam >> $out/varReads/bam/bams_$prefix.txt";}
	
	system "samtools merge -b $out/varReads/bam/bams_$prefix.txt $out/varReads/bam/var_$prefix.bam";
	
	system "samtools fastq -1 $out/varReads/1/$prefix.1.fq -2 $out/varReads/2/$prefix.2.fq -0 /dev/null -s $out/varReads/singleton/$prefix.fq -n $out/varReads/bam/var_$prefix.bam";		
    }
}
