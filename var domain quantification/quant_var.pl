#!/usr/bin/perl
##
##
##use warnings;
use strict;

if ($#ARGV != 2) {
    print "usage: Enter full path to target genome directory, output directory, directory with reads \ne.g. countRNASeqReads.pl /home/genome /home/out /home/reads/1 /home/reads/2 /home/file.gtf";
    exit;
}
my $GenomeDir = $ARGV[0];
my $out = $ARGV[1];
my $pair1 = $ARGV[2];

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
my @Readfiles1 = ();
@Readfiles1 = getReads("$pair1");

align($GenomeDir, "$pair1", "$out", \@Readfiles1);
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
        unless (-d "$topOut/quant")
        {system "mkdir $topOut/quant";}
        print "\tSalmon quant files stored in $topOut/quant\n";
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
###########################################################
sub align
{
    my $GenomeDir = shift;
    my $read1Dir = shift;
    my $out = shift;
    my $pair1Reads = shift;
    my @Reads1 = @{$pair1Reads};
    my $i;
    my $size = @Reads1;
    my @PRE;
    my $prefix;

    for ($i = 0; $i < $size; $i++)
    {
        @PRE = split(/\./, $Reads1[$i]);
        $prefix = $PRE[0];

        print "\n$prefix\n";
        system "mkdir $out/quant/$prefix";

        my $read1 = "$read1Dir/$prefix" . ".fq";
		
	system "salmon quant -i $GenomeDir -l A -r $read1 -p 8 --validateMappings -o $out/quant/$prefix";
    }
}

