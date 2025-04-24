#!/usr/bin/perl -w

###################################################################################################
###############parse the sequences (-200/+50) based on the polyA site information##################
###################################################################################################
my $gmapdir="/usr/local/genome/gmap/gmap-2013-05-09/bin/";
my $database='t3702.tair10.gmap20110831.k12';                ###type your database name here;
my $order='';
my $strand='';
my $chr='';
my $cmd='';
my $site=0;
my $start=0;
my $end=0;
my $cigar='';

my $contents1='';
my $contents2='';
my $contents3='';
my $contents4='';
my $contents5='';
my $contents6='';
my $contents='';


my $filename=shift;
open(F,"<$filename")||die "$!";
open(OF, ">$filename.clean")||die "$!";
while(my $line=<F>)
{
	chomp($line);
    my @temp=split(/\t/, $line);
	if($temp[0] eq '@HD')
	{
		next;
	}
    if($temp[0] eq '@PG')
	{
		next;
	}
	$order=$temp[0];
	$strand=$temp[1];
    $chr=$temp[2];
    $site=$temp[3];
    $contents1=$temp[4];
    $cigar=$temp[5];
    $contents2=$temp[6];
    $contents3=$temp[7];
    $contents4=$temp[8];
    $seq=$temp[9];
    $contents5=$temp[10];
    if($strand eq '16')
    {                            
		 my @cigar_temp = split(//, $cigar);
		 my $insert=0;
		 my $deletion=0;

		 for(my $i=0;$i<@cigar_temp;$i++)
		 {
		 	
		 	if($cigar_temp[$i] eq 'I')
		 	{
		 		$insert=$insert + int($cigar_temp[$i-1]);
		 	}
		 	if($cigar_temp[$i] eq 'D')
		 	{
		 		$deletion=$deletion+int($cigar_temp[$i-1]);
		 	}
		 }
		 
		 $site=$site + length($seq) -$insert +$deletion; 
    }
	my $internal_primer_count=0;
	my $internal_primer_count_rc=0;
    if($strand eq '0')
    {                             ##positive strand;
		$start=$site-15;
		$end=$site;
		$cmd = $gmapdir."get-genome -d $database $chr:$start..$end";       
		$contents=`$cmd`;
    	$contents=~s/>.+\d+//g; # removing header info in contents
    	$contents=~s/[\n\r\f\^M]//g; #removing unwanted char
    	
    	my @internal_primer = split(//, $contents);
    	
    	for(my $j=0;$j<@internal_primer;$j++)
		{
			
			if($internal_primer[$j] eq 'T')
			{
			 	$internal_primer_count++;
			}
			
		}
    	
    }
    else
    {                                          ##negative strand;
		$start=$site;
		$end=$site+15;
		$cmd = $gmapdir."get-genome -d $database $chr:$start..$end";  
		$contents=`$cmd`;
    	$contents=~s/>.+\d+//g; # removing header info in contents
    	$contents=~s/[\n\r\f\^M]//g; #removing unwanted char
    	
    	my @internal_primer = split(//, $contents);
    	
    	for(my $j=0;$j<@internal_primer;$j++)
		{
			
			if($internal_primer[$j] eq 'A')
			{
			 	$internal_primer_count++;
			}
			
		}
		
    }

    if($internal_primer_count<9)
    {
    	print OF $order, "\t", $strand, "\t", $chr, "\t", $site, "\t", $contents1, "\t", $cigar, "\t", $contents2, "\t", $contents3, "\t", $contents4, "\t", $seq, "\t", $contents5, "\n";
    }

}
close(OF);
close(F);