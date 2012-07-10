#!/usr/bin/perl
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: Sat 07 Jul 2012 09:20:30 PM CDT
#Vanderbilt Center for Quantitative Sciences
#############################################
use strict;
use warnings;
use Getopt::Long;
#############################################
my (
	$inputVcfFile,       $dbSNPVcfFile,            $nonsynonymousFile,
	$dbSNPVcfFileFormat, $nonsynonymousFileFormat, $nsoft,
	$help,               $qualFilter,              $filterFilter,
	$gqFilter,           $dpFilter
);

my $output = ".";

my $commandline = join " ", ( "perl vcfStat.pl", @ARGV );

#print $commandline,"\n\n"; exit(1);
exit(1)
  unless GetOptions(
	"i|input=s"         => \$inputVcfFile,
	"d|dbsnp=s"         => \$dbSNPVcfFile,
	"df|dbsnpformat"    => \$dbSNPVcfFileFormat,
	"n|nonsynonymous=s" => \$nonsynonymousFile,
	"nf"                => \$nonsynonymousFileFormat,
	"h|help"            => \$help,
	"qual=f"            => \$qualFilter,
	"filter=s"          => \$filterFilter,
	"gq=f"              => \$gqFilter,
	"dp=f"              => \$dpFilter,
	"o|out=s"           => \$output
  );

################################################
my $usage = <<USAGE;
*************************************************************************************************
Description: To get a summary of SNPs from a VCF file. Information includes Ti/Tv, novel SNPs etc.

Usage: perl vcfStat.pl -i/--input <inputfile> -d/--dbsnp <dbsnp file> -n/--nonsynonymous <non-synonymousfile> -h/--help -nf --qual --filter --gq --dp
          -i/--input          input vcffile (*required)
          -d/--dbsnp          input dbsnp file either in vcf format or a 'list'(chr pos [rsid] [ref] [alt]), default is vcf format (*required)
          -df/--dbsnpformat   if this option is provided, it means the dbsnp file in a 'list' format (optional)
          -n/--nonsynoymous   non-synonymous snps file annotated by ANNOVAR. The default format assumes the 2nd, 4th and 5th columns
                              are 'snptype(synonymous, non-synonymous etc.)','chr' and 'pos'. If option provided with '-nf', it means 
                              the input is in a 'list' format (chr pos [rsid] [ref] [alt]) (*required)
          -nf                 if this option is provided, it means the non-synonymous file in a 'list' format (optional)
          -nsoft              if this option provided, 'non-synonymous' SNPs not include those 'stopgain' and 'stoploss' ones. Default inclues those two types. (optional)
          --qual              'qual' filter for the input vcf file (qual < given cutoff will be skipped). The qual filter is the 6th column in a vcf.(optional)
          --filter            'filter' filter for the input vcf file, only FILTER name eq to 'filter' in the 7th columns will be kept. Default is 'PASS'. 
                              Give '*' to make no filter. Multile filter names could be give by 'filter1|filter2|filter3'.(optional)
          --gq                Genotype quality filter for each sample. If a given samples' GQ less than cutoff, the genotype will be masked as './.'. (optional)
          --dp                Depth filter for each sample. If a given samples' DP less than cutoff, the genotype will be masked as './.'. (optional)
          -o/--out            output filename for the report, default is vcf_report.html (optional)
          -h/--help           help page
USAGE

########################################
#Check
if ($help) {
	print $usage;
	exit(0);
}

# 1) input
if ( !defined($inputVcfFile) ) {
	print "[Error] input vcf file is not provided yet\n";
	print $usage;
	exit(1);
}
elsif ( !-e $inputVcfFile ) {
	print "[Error] input vcf file '$inputVcfFile' is not exists\n";
	print $usage;
	exit(1);
}

# 2 dbsnp
if ( !defined($dbSNPVcfFile) ) {
	print "[Error] input dbSNP file is not provided yet\n";
	print $usage;
	exit(1);
}
elsif ( !-e $dbSNPVcfFile ) {
	print "[Error] input dbSNP file '$dbSNPVcfFile' is not exists\n";
	print $usage;
	exit(1);
}

# 3 non-synonymous
if ( !defined($nonsynonymousFile) ) {
	print "[Error] input non-synonymous file is not provided yet\n";
	print $usage;
	exit(1);
}
elsif ( !-e $nonsynonymousFile ) {
	print
	  "[Error] input non-synonymous file '$nonsynonymousFile' is not exists\n";
	print $usage;
	exit(1);
}

#end check
####################################

my %dbsnp;            #keys are chr\tpos
my %nonsynonymous;    #keys are chr\tpos

#######################
#Main

#1 load dbsnp
info("Loading dbSNP information ...");
loadDBSNP( $dbSNPVcfFile, \%dbsnp, $dbSNPVcfFileFormat );
info( "There are " . scalar( keys %dbsnp ) . " SNPs in dbSNP" );

# 2 load non-synonymous
info("Loading non-synonymous SNPs' information ...");
loadNonSynonymous( $nonsynonymousFile, \%nonsynonymous,
	$nonsynonymousFileFormat, $nsoft );
info( "There are " . scalar( keys %nonsynonymous ) . " non-synonymous SNPs" );

# 3 parse vcf
info("Loading input vcf file ...");
my %data;    #To store data
my ( $chrref, $sampleref ) = parseVCF( $inputVcfFile, \%data );

# 4 calculate titv in the $ref, and initialize some value
normData( \%data, $chrref, $sampleref );

info("Write out result into data sub-folder...");
writeData( \%data );

info("Generate report ...\n");
writeHtml( \%data );

#Load dbSNP into hash
sub loadDBSNP {
	my ( $in, $ref, $flag ) = @_;
	open( IN, $in ) or die $!;
	if ($flag) {    #means a 'list' format file
		while (<IN>) {
			s/\r|\n//g;
			next if (/^#|^$/);
			my ( $chr, $pos, $rsid, $refallele, $altallele ) = split /\s+/;
			if ( $pos !~ /^\d/ ) {
				print "[Error] Input dbSNP file is not in 'list' format \n";
				print $usage;
				exit(1);
			}
			my $key = join "\t", ( $chr, $pos );
			$ref->{$key} = undef;
		}
	}
	else {
		while (<IN>) {
			s/\r|\n//g;
			next if (/^#|^$/);
			my (
				$chr,  $pos,    $rsid, $refallele, $altallele,
				$qual, $filter, $info, $format,    $details
			) = split /\s+/;
			if ( $pos !~ /^\d/ ) {
				print "[Error] Input dbSNP file is not in 'vcf' format \n";
				print $usage;
				exit(1);
			}
			my $key = join "\t", ( $chr, $pos );
			$ref->{$key} = undef;
		}
	}
	close IN;
}

#Load non-synonymous SNPs
sub loadNonSynonymous {
	my ( $in, $ref, $flag, $flag2 ) = @_;
	open( IN, $in ) or die $!;
	if ($flag) {    #means a 'list' format file
		while (<IN>) {
			s/\r|\n//g;
			next if (/^#|^$/);
			my ( $chr, $pos, $rsid, $refallel, $altallele ) = split /\s+/;
			if ( $pos !~ /^\d/ ) {
				print
"[Error] Input non-synonymous file is not in 'list' format \n";
				print $usage;
				exit(1);
			}
			my $key = join "\t", ( $chr, $pos );
			$ref->{$key} = undef;
		}
	}
	else {
		while (<IN>) {
			s/\r|\n//g;
			next if (/^#|^$/);
			my ( $line, $type, $details, $chr, $pos, @others ) = split /\s+/;
			if ( $pos !~ /^\d/ ) {
				print
				  "[Error] Input non-synonymous file is not in right format \n";
				print $usage;
				exit(1);
			}

			my $key = join "\t", ( $chr, $pos );
			if ($flag2) {    #only snp masked as nonsynonymous SNV
				if ( $type eq 'nonsynonymous SNV' ) {
					$ref->{$key} = undef;
				}
			}
			else {
				if (   $type eq 'nonsynonymous SNV'
					|| $type eq 'stopgain SNV'
					|| $type eq 'stopgain SNV' )
				{
					$ref->{$key} = undef;
				}
			}
		}
	}
	close IN;
}

#$ref structure

=ref structure
ref={'basic'=>{
				'overall'=>{'total'=>number,'remained'=>number},
				'snp'=>{'total'=>number,'remained'=>number},
				'others'=>{'total'=>number,'remained'=>number}, #may be indel
			 },
			 
	 'snp'=>{
				'allchr'=>{'overall'=>{'number'=>number,'ti'=>number,'tv'=>number},
						   'novel'=>{'number'=>number,'ti'=>number,'tv'=>number},
						   'non-synonymous'=>{'number'=>number,'ti'=>number,'tv'=>number}
						   },
				'perchr'=>{
							'chr1'=>{'overall'=>{'number'=>number,'ti'=>number,'tv'=>number},
									   'novel'=>{'number'=>number,'ti'=>number,'tv'=>number},
									   'non-synonymous'=>{'number'=>number,'ti'=>number,'tv'=>number}
									},
							'chr2'=>{'overall'=>{'number'=>number,'ti'=>number,'tv'=>number},
									   'novel'=>{'number'=>number,'ti'=>number,'tv'=>number},
									   'non-synonymous'=>{'number'=>number,'ti'=>number,'tv'=>number}
									}
						   }
			}
	'samples'=>{
				'sample'=>{
						'allchr'=>{
						   'overall'=>{'number'=>number,'ti'=>number,'tv'=>number},
						   'novel'=>{'number'=>number,'ti'=>number,'tv'=>number},
						   'non-synonymous'=>{'number'=>number,'ti'=>number,'tv'=>number},
							},
						'perchr'=>{
							'chr1'=>{'overall'=>{'number'=>number,'ti'=>number,'tv'=>number},
									   'novel'=>{'number'=>number,'ti'=>number,'tv'=>number},
									   'non-synonymous'=>{'number'=>number,'ti'=>number,'tv'=>number}
									},
							'chr2'=>{'overall'=>{'number'=>number,'ti'=>number,'tv'=>number},
									   'novel'=>{'number'=>number,'ti'=>number,'tv'=>number},
									   'non-synonymous'=>{'number'=>number,'ti'=>number,'tv'=>number}
									}
						}
				}
			}
}
=cut

sub parseVCF {
	my ( $in, $ref ) = @_;
	open( IN, $in ) or die $!;
	my %header;    #key is sample name, value is its index, 0-based
	my %chr;
	while (<IN>) {
		s/\r|\n//g;
		next if (/^##|^$/);
		if (/^#/) {
			my (
				$chr,  $pos,    $id,   $refAllele, $altAllele,
				$qual, $filter, $info, $format,    @samples
			) = split "\t";
			for ( my $i = 0 ; $i < @samples ; $i++ ) {
				$header{ $samples[$i] } = $i;
			}
			next;
		}

		#Parse real data
		my (
			$chr,  $pos,    $id,   $refAllele, $altAllele,
			$qual, $filter, $info, $format,    @genotypes
		) = split "\t";
		$chr{$chr}++;
		$ref->{'basic'}->{'overall'}->{'total'}++;
		my $issnp = isSNP( $refAllele, $altAllele );
		if ($issnp) {
			$ref->{'basic'}->{'snp'}->{'total'}++;
		}
		else {
			$ref->{'basic'}->{'others'}->{'total'}++;
		}

#If this postion pass the filters user defined, then continue, otherwise go to next position
		next if ( $qualFilter && $qual < $qualFilter );

		if ($filterFilter) {
			my @tmp = split /\|/, $filterFilter;
			my $flag = 0;    #default not passed
			foreach my $t (@tmp) {
				if ( uc($t) eq uc($filter) ) {
					$flag = 1;    #passed the filter
					next;
				}
			}
			next unless ($flag);    #flag=0, not pass the filter
		}

		my @formats = split ":", $format;
		my %temp;
		map { $temp{$_} = undef } @formats;
		my @newGenotypes;           #Only contains genotype value
		my %newGenotypes;
		if (   ( $gqFilter && exists( $temp{'GQ'} ) )
			|| ( $dpFilter && exists( $temp{'DP'} ) ) )
		{                           #user definded one of these two filters
			foreach my $geno (@genotypes) {
				my %map;
				@map{@formats} = split ":", $geno;
				my $gt = $map{'GT'};
				if ( $gt =~ /\./ ) {
					$gt = "./.";
				}
				else {
					if (   $gqFilter
						&& exists( $map{'GQ'} )
						&& $map{'GQ'} < $gqFilter )
					{
						$gt = "./.";
					}
					if (   $dpFilter
						&& exists( $map{'DP'} )
						&& $map{'DP'} < $dpFilter )
					{
						$gt = "./.";
					}
				}
				$gt =~ s/\|/\//g;
				push @newGenotypes, $gt;
				$newGenotypes{$gt}++;
			}
		}
		else {
			foreach my $geno (@genotypes) {
				my %map;
				@map{@formats} = split ":", $geno;
				my $gt = $map{'GT'};
				if ( $gt =~ /\./ ) {
					$gt = "./.";
				}
				$gt =~ s/\|/\//g;
				push @newGenotypes, $gt;
				$newGenotypes{$gt}++;
			}
		}

#If across all the samples, there is no alternative allele called, then skip this position
		my $missingGenotype = 0;
		my $refHomoGenotype = 0;
		if ( exists( $newGenotypes{'./.'} ) ) {
			$missingGenotype = $newGenotypes{'./.'};
		}
		if ( exists( $newGenotypes{'0/0'} ) ) {
			$refHomoGenotype = $newGenotypes{'0/0'};
		}
		next
		  if ( $missingGenotype + $refHomoGenotype == scalar(@newGenotypes) );

		#remained
		my $key = join "\t", ( $chr, $pos );

		my $isti = isTi( $refAllele, $altAllele );
		$ref->{'basic'}->{'overall'}->{'remained'}++;
		unless ($issnp) {
			$ref->{'basic'}->{'others'}->{'remained'}++;
		}
		else {
			$ref->{'basic'}->{'snp'}->{'remained'}++;

		#Do other more things;
		#total snp/ti snp/ tv snp/ numbers on all chromosome and each chromosome
			$ref->{'snp'}->{'allchr'}->{'overall'}->{'number'}++;
			$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'number'}++;
			if ($isti) {
				$ref->{'snp'}->{'allchr'}->{'overall'}->{'ti'}++;
				$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'ti'}++;
			}
			else {
				$ref->{'snp'}->{'allchr'}->{'overall'}->{'tv'}++;
				$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'tv'}++;
			}

			if ( $dbSNPVcfFile && %dbsnp && !exists( $dbsnp{$key} ) )
			{    #novel snp
				$ref->{'snp'}->{'allchr'}->{'novel'}->{'number'}++;
				$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'number'}++;
				if ($isti) {
					$ref->{'snp'}->{'allchr'}->{'novel'}->{'ti'}++;
					$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'ti'}++;
				}
				else {
					$ref->{'snp'}->{'allchr'}->{'novel'}->{'tv'}++;
					$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'tv'}++;
				}
			}
			if (   $nonsynonymousFile
				&& %nonsynonymous
				&& exists( $nonsynonymous{$key} ) )
			{
				$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'number'}++;
				$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}
				  ->{'number'}++;
				if ($isti) {
					$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'ti'}++;
					$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}
					  ->{'ti'}++;
				}
				else {
					$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'tv'}++;
					$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}
					  ->{'tv'}++;
				}
			}

			if (   $nonsynonymousFile
				&& %nonsynonymous
				&& exists( $nonsynonymous{$key} )
				&& !exists( $dbsnp{$key} ) )
			{
				$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}
				  ->{'number'}++;
				$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
				  ->{'number'}++;
				if ($isti) {
					$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}
					  ->{'ti'}++;
					$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
					  ->{'ti'}++;
				}
				else {
					$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}
					  ->{'tv'}++;
					$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
					  ->{'tv'}++;
				}
			}
		}

		#loop each sample
		foreach my $sample ( sort keys %header ) {
			my $index    = $header{$sample};
			my $genotype = $newGenotypes[$index];
			next
			  if ( $genotype !~ /1/ )
			  ; #only genotype contains alternatvie allele is considered as a variants

			$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'number'}++;
			$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
			  ->{'number'}++;
			if ($isti) {
				$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'ti'}++;
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
				  ->{'ti'}++;
			}
			else {
				$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'tv'}++;
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
				  ->{'tv'}++;
			}

			if ( $dbSNPVcfFile && %dbsnp && !exists( $dbsnp{$key} ) )
			{    #novel snp
				$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}
				  ->{'number'}++;
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
				  ->{'number'}++;
				if ($isti) {
					$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}
					  ->{'ti'}++;
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
					  ->{'ti'}++;
				}
				else {
					$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}
					  ->{'tv'}++;
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
					  ->{'tv'}++;
				}
			}

			if (   $nonsynonymousFile
				&& %nonsynonymous
				&& exists( $nonsynonymous{$key} ) )
			{
				$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
				  ->{'number'}++;
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'non-synonymous'}->{'number'}++;
				if ($isti) {
					$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
					  ->{'ti'}++;
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'non-synonymous'}->{'ti'}++;
				}
				else {
					$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
					  ->{'tv'}++;
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'non-synonymous'}->{'tv'}++;
				}
			}

			if (   $nonsynonymousFile
				&& %nonsynonymous
				&& exists( $nonsynonymous{$key} )
				&& !exists( $dbsnp{$key} ) )
			{
				$ref->{'samples'}->{$sample}->{'allchr'}
				  ->{'novel_non-synonymous'}->{'number'}++;
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'novel_non-synonymous'}->{'number'}++;
				if ($isti) {
					$ref->{'samples'}->{$sample}->{'allchr'}
					  ->{'novel_non-synonymous'}->{'ti'}++;
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'novel_non-synonymous'}->{'ti'}++;
				}
				else {
					$ref->{'samples'}->{$sample}->{'allchr'}
					  ->{'novel_non-synonymous'}->{'tv'}++;
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'novel_non-synonymous'}->{'tv'}++;
				}
			}
		}
	}
	close IN;
	return ( \%chr, \%header );
}

#calculate titv in the $ref, and initialize some value
sub normData {
	my ( $ref, $chrref, $sampleref ) = @_;

	#1) basic
	$ref->{'basic'}->{'overall'}->{'total'} = 0
	  if ( !exists( $ref->{'basic'}->{'overall'}->{'total'} ) );
	$ref->{'basic'}->{'overall'}->{'remained'} = 0
	  if ( !exists( $ref->{'basic'}->{'overall'}->{'remained'} ) );
	$ref->{'basic'}->{'snp'}->{'total'} = 0
	  if ( !exists( $ref->{'basic'}->{'snp'}->{'total'} ) );
	$ref->{'basic'}->{'snp'}->{'remained'} = 0
	  if ( !exists( $ref->{'basic'}->{'snp'}->{'remained'} ) );
	$ref->{'basic'}->{'others'}->{'total'} = 0
	  if ( !exists( $ref->{'basic'}->{'others'}->{'total'} ) );
	$ref->{'basic'}->{'others'}->{'remained'} = 0
	  if ( !exists( $ref->{'basic'}->{'others'}->{'remained'} ) );

	$ref->{'basic'}->{'overall'}->{'filtered'} =
	  $ref->{'basic'}->{'overall'}->{'total'} -
	  $ref->{'basic'}->{'overall'}->{'remained'};
	$ref->{'basic'}->{'snp'}->{'filtered'} =
	  $ref->{'basic'}->{'snp'}->{'total'} -
	  $ref->{'basic'}->{'overall'}->{'remained'};
	$ref->{'basic'}->{'others'}->{'filtered'} =
	  $ref->{'basic'}->{'others'}->{'total'} -
	  $ref->{'basic'}->{'others'}->{'remained'};

	#2) snp information
	$ref->{'snp'}->{'allchr'}->{'overall'}->{'number'} = 0
	  if ( !exists( $ref->{'snp'}->{'allchr'}->{'overall'}->{'number'} ) );
	$ref->{'snp'}->{'allchr'}->{'overall'}->{'ti'} = 0
	  if ( !exists( $ref->{'snp'}->{'allchr'}->{'overall'}->{'ti'} ) );
	$ref->{'snp'}->{'allchr'}->{'overall'}->{'tv'} = 0
	  if ( !exists( $ref->{'snp'}->{'allchr'}->{'overall'}->{'tv'} ) );
	if ( $ref->{'snp'}->{'allchr'}->{'overall'}->{'tv'} == 0 ) {
		$ref->{'snp'}->{'allchr'}->{'overall'}->{'titv'} = 'NA';
	}
	else {
		$ref->{'snp'}->{'allchr'}->{'overall'}->{'titv'} =
		  $ref->{'snp'}->{'allchr'}->{'overall'}->{'ti'} /
		  $ref->{'snp'}->{'allchr'}->{'overall'}->{'tv'};
	}

	if ( $dbSNPVcfFile && %dbsnp ) {
		$ref->{'snp'}->{'allchr'}->{'novel'}->{'number'} = 0
		  if ( !exists( $ref->{'snp'}->{'allchr'}->{'novel'}->{'number'} ) );
		$ref->{'snp'}->{'allchr'}->{'novel'}->{'ti'} = 0
		  if ( !exists( $ref->{'snp'}->{'allchr'}->{'novel'}->{'ti'} ) );
		$ref->{'snp'}->{'allchr'}->{'novel'}->{'tv'} = 0
		  if ( !exists( $ref->{'snp'}->{'allchr'}->{'novel'}->{'tv'} ) );
		if ( $ref->{'snp'}->{'allchr'}->{'novel'}->{'tv'} == 0 ) {
			$ref->{'snp'}->{'allchr'}->{'novel'}->{'titv'} = 'NA';
		}
		else {
			$ref->{'snp'}->{'allchr'}->{'novel'}->{'titv'} =
			  $ref->{'snp'}->{'allchr'}->{'novel'}->{'ti'} /
			  $ref->{'snp'}->{'allchr'}->{'novel'}->{'tv'};
		}
	}
	else {
		$ref->{'snp'}->{'allchr'}->{'novel'}->{'number'} = 'NA';
		$ref->{'snp'}->{'allchr'}->{'novel'}->{'titv'}   = 'NA';
	}
	if ( $nonsynonymousFile && %nonsynonymous ) {
		$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'number'} = 0
		  if (
			!exists(
				$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'number'}
			)
		  );
		$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'ti'} = 0
		  if (
			!exists( $ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'ti'} ) );
		$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'tv'} = 0
		  if (
			!exists( $ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'tv'} ) );
		if ( $ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'tv'} == 0 ) {
			$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'titv'} = 'NA';
		}
		else {
			$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'titv'} =
			  $ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'ti'} /
			  $ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'tv'};
		}
	}
	else {
		$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'number'} = 'NA';
		$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'titv'}   = 'NA';
	}

	if ( $nonsynonymousFile && %nonsynonymous && %dbsnp ) {
		$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'number'} = 0
		  if (
			!exists(
				$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'number'}
			)
		  );
		$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'ti'} = 0
		  if (
			!exists(
				$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'ti'}
			)
		  );
		$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'tv'} = 0
		  if (
			!exists(
				$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'tv'}
			)
		  );
		if ( $ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'tv'} == 0 )
		{
			$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'titv'} =
			  'NA';
		}
		else {
			$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'titv'} =
			  $ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'ti'} /
			  $ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'tv'};
		}
	}
	else {
		$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'number'} = 'NA';
		$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'titv'}   = 'NA';
	}

	#for each chromosome
	foreach my $chr ( keys %{$chrref} ) {
		$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'number'} = 0
		  if (
			!exists(
				$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'number'}
			)
		  );
		$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'ti'} = 0
		  if (
			!exists( $ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'ti'} ) );
		$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'tv'} = 0
		  if (
			!exists( $ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'tv'} ) );
		if ( $ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'tv'} == 0 ) {
			$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'titv'} = 'NA';
		}
		else {
			$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'titv'} =
			  $ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'ti'} /
			  $ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'tv'};
		}

		if ( $dbSNPVcfFile && %dbsnp ) {
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'number'} = 0
			  if (
				!exists(
					$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'number'}
				)
			  );
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'ti'} = 0
			  if (
				!exists( $ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'ti'} )
			  );
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'tv'} = 0
			  if (
				!exists( $ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'tv'} )
			  );
			if ( $ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'tv'} == 0 ) {
				$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'titv'} = 'NA';
			}
			else {
				$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'titv'} =
				  $ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'ti'} /
				  $ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'tv'};
			}
		}
		else {
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'number'} = 'NA';
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'titv'}   = 'NA';
		}
		if ( $nonsynonymousFile && %nonsynonymous ) {
			$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}->{'number'} =
			  0
			  if (
				!exists(
					$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}
					  ->{'number'}
				)
			  );
			$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}->{'ti'} = 0
			  if (
				!exists(
					$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}
					  ->{'ti'}
				)
			  );
			$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}->{'tv'} = 0
			  if (
				!exists(
					$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}
					  ->{'tv'}
				)
			  );
			if (
				$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}->{'tv'} ==
				0 )
			{
				$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}
				  ->{'titv'} = 'NA';
			}
			else {
				$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}
				  ->{'titv'} =
				  $ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}
				  ->{'ti'} /
				  $ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}->{'tv'};
			}
		}
		else {
			$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}->{'number'} =
			  'NA';
			$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}->{'titv'} =
			  'NA';
		}

		if ( $nonsynonymousFile && %nonsynonymous && %dbsnp ) {
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
			  ->{'number'} = 0
			  if (
				!exists(
					$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
					  ->{'number'}
				)
			  );
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
			  ->{'ti'} = 0
			  if (
				!exists(
					$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
					  ->{'ti'}
				)
			  );
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
			  ->{'tv'} = 0
			  if (
				!exists(
					$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
					  ->{'tv'}
				)
			  );
			if ( $ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
				->{'tv'} == 0 )
			{
				$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
				  ->{'titv'} = 'NA';
			}
			else {
				$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
				  ->{'titv'} =
				  $ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
				  ->{'ti'} /
				  $ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
				  ->{'tv'};
			}
		}
		else {
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
			  ->{'number'} = 'NA';
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
			  ->{'titv'} = 'NA';
		}
	}

	#3) sample level
	foreach my $sample ( keys %{$sampleref} ) {
		$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'number'} = 0
		  if (
			!exists(
				$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}
				  ->{'number'}
			)
		  );
		$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'ti'} = 0
		  if (
			!exists(
				$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'ti'}
			)
		  );
		$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'tv'} = 0
		  if (
			!exists(
				$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'tv'}
			)
		  );
		if (
			$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'tv'} == 0 )
		{
			$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'titv'} =
			  'NA';
		}
		else {
			$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'titv'} =
			  $ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'ti'} /
			  $ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'tv'};
		}

		if ( $dbSNPVcfFile && %dbsnp ) {
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'number'} = 0
			  if (
				!exists(
					$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}
					  ->{'number'}
				)
			  );
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'ti'} = 0
			  if (
				!exists(
					$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'ti'}
				)
			  );
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'tv'} = 0
			  if (
				!exists(
					$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'tv'}
				)
			  );
			if ( $ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'tv'} ==
				0 )
			{
				$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'titv'} =
				  'NA';
			}
			else {
				$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'titv'} =
				  $ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'ti'} /
				  $ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'tv'};
			}
		}
		else {
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'number'} =
			  'NA';
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'titv'} =
			  'NA';
		}
		if ( $nonsynonymousFile && %nonsynonymous ) {
			$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
			  ->{'number'} = 0
			  if (
				!exists(
					$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
					  ->{'number'}
				)
			  );
			$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
			  ->{'ti'} = 0
			  if (
				!exists(
					$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
					  ->{'ti'}
				)
			  );
			$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
			  ->{'tv'} = 0
			  if (
				!exists(
					$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
					  ->{'tv'}
				)
			  );
			if ( $ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
				->{'tv'} == 0 )
			{
				$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
				  ->{'titv'} = 'NA';
			}
			else {
				$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
				  ->{'titv'} =
				  $ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
				  ->{'ti'} /
				  $ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
				  ->{'tv'};
			}
		}
		else {
			$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
			  ->{'number'} = 'NA';
			$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
			  ->{'titv'} = 'NA';
		}

		if ( $nonsynonymousFile && %nonsynonymous && %dbsnp ) {
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel_non-synonymous'}
			  ->{'number'} = 0
			  if (
				!exists(
					$ref->{'samples'}->{$sample}->{'allchr'}
					  ->{'novel_non-synonymous'}->{'number'}
				)
			  );
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel_non-synonymous'}
			  ->{'ti'} = 0
			  if (
				!exists(
					$ref->{'samples'}->{$sample}->{'allchr'}
					  ->{'novel_non-synonymous'}->{'ti'}
				)
			  );
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel_non-synonymous'}
			  ->{'tv'} = 0
			  if (
				!exists(
					$ref->{'samples'}->{$sample}->{'allchr'}
					  ->{'novel_non-synonymous'}->{'tv'}
				)
			  );
			if ( $ref->{'samples'}->{$sample}->{'allchr'}
				->{'novel_non-synonymous'}->{'tv'} == 0 )
			{
				$ref->{'samples'}->{$sample}->{'allchr'}
				  ->{'novel_non-synonymous'}->{'titv'} = 'NA';
			}
			else {
				$ref->{'samples'}->{$sample}->{'allchr'}
				  ->{'novel_non-synonymous'}->{'titv'} =
				  $ref->{'samples'}->{$sample}->{'allchr'}
				  ->{'novel_non-synonymous'}->{'ti'} /
				  $ref->{'samples'}->{$sample}->{'allchr'}
				  ->{'novel_non-synonymous'}->{'tv'};
			}
		}
		else {
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel_non-synonymous'}
			  ->{'number'} = 'NA';
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel_non-synonymous'}
			  ->{'titv'} = 'NA';
		}

		#loop each chromosome
		foreach my $chr ( keys %{$chrref} ) {
			$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
			  ->{'number'} = 0
			  if (
				!exists(
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'overall'}->{'number'}
				)
			  );
			$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
			  ->{'ti'} = 0
			  if (
				!exists(
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'overall'}->{'ti'}
				)
			  );
			$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
			  ->{'tv'} = 0
			  if (
				!exists(
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'overall'}->{'tv'}
				)
			  );
			if ( $ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
				->{'tv'} == 0 )
			{
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
				  ->{'titv'} = 'NA';
			}
			else {
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
				  ->{'titv'} =
				  $ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
				  ->{'ti'} /
				  $ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
				  ->{'tv'};
			}

			if ( $dbSNPVcfFile && %dbsnp ) {
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
				  ->{'number'} = 0
				  if (
					!exists(
						$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
						  ->{'novel'}->{'number'}
					)
				  );
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
				  ->{'ti'} = 0
				  if (
					!exists(
						$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
						  ->{'novel'}->{'ti'}
					)
				  );
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
				  ->{'tv'} = 0
				  if (
					!exists(
						$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
						  ->{'novel'}->{'tv'}
					)
				  );
				if ( $ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
					->{'tv'} == 0 )
				{
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
					  ->{'titv'} = 'NA';
				}
				else {
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
					  ->{'titv'} =
					  $ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'novel'}->{'ti'} /
					  $ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'novel'}->{'tv'};
				}
			}
			else {
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
				  ->{'number'} = 'NA';
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
				  ->{'titv'} = 'NA';
			}
			if ( $nonsynonymousFile && %nonsynonymous ) {
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'non-synonymous'}->{'number'} = 0
				  if (
					!exists(
						$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
						  ->{'non-synonymous'}->{'number'}
					)
				  );
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'non-synonymous'}->{'ti'} = 0
				  if (
					!exists(
						$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
						  ->{'non-synonymous'}->{'ti'}
					)
				  );
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'non-synonymous'}->{'tv'} = 0
				  if (
					!exists(
						$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
						  ->{'non-synonymous'}->{'tv'}
					)
				  );
				if ( $ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					->{'non-synonymous'}->{'tv'} == 0 )
				{
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'non-synonymous'}->{'titv'} = 'NA';
				}
				else {
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'non-synonymous'}->{'titv'} =
					  $ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'non-synonymous'}->{'ti'} /
					  $ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'non-synonymous'}->{'tv'};
				}
			}
			else {
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'non-synonymous'}->{'number'} = 'NA';
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'non-synonymous'}->{'titv'} = 'NA';
			}

			if ( $nonsynonymousFile && %nonsynonymous && %dbsnp ) {
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'novel_non-synonymous'}->{'number'} = 0
				  if (
					!exists(
						$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
						  ->{'novel_non-synonymous'}->{'number'}
					)
				  );
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'novel_non-synonymous'}->{'ti'} = 0
				  if (
					!exists(
						$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
						  ->{'novel_non-synonymous'}->{'ti'}
					)
				  );
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'novel_non-synonymous'}->{'tv'} = 0
				  if (
					!exists(
						$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
						  ->{'novel_non-synonymous'}->{'tv'}
					)
				  );
				if ( $ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					->{'novel_non-synonymous'}->{'tv'} == 0 )
				{
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'novel_non-synonymous'}->{'titv'} = 'NA';
				}
				else {
					$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'novel_non-synonymous'}->{'titv'} =
					  $ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'novel_non-synonymous'}->{'ti'} /
					  $ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
					  ->{'novel_non-synonymous'}->{'tv'};
				}
			}
			else {
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'novel_non-synonymous'}->{'number'} = 'NA';
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'novel_non-synonymous'}->{'titv'} = 'NA';
			}
		}
	}
}

sub writeData {
	my ($ref) = @_;
	if ( !-d "data" ) {
		mkdir "data";
	}

	#1) Write out basic information
	open( OUT, ">data/basic_info.txt" ) or die $!;
	print OUT "#Types\tTotalVariants\tFilteredOut\n";
	print OUT join "\t",
	  (
		"Overalll", $ref->{'basic'}->{'overall'}->{'total'},
		$ref->{'basic'}->{'overall'}->{'total'} -
		  $ref->{'basic'}->{'overall'}->{'remained'}
	  );
	print OUT "\n";
	print OUT join "\t",
	  (
		"SNPs", $ref->{'basic'}->{'snp'}->{'total'},
		$ref->{'basic'}->{'snp'}->{'total'} -
		  $ref->{'basic'}->{'snp'}->{'remained'}
	  );
	print OUT "\n";
	print OUT join "\t",
	  (
		"Others", $ref->{'basic'}->{'others'}->{'total'},
		$ref->{'basic'}->{'others'}->{'total'} -
		  $ref->{'basic'}->{'others'}->{'remained'}
	  );
	print OUT "\n";
	close OUT;

	#2) Write Out SNP information
	open( OUT, ">data/snp_info.txt" ) or die $!;
	print OUT
">#Types\tOverallSNPs\tOverallTiTv\tNovelSNPs\tNovelTiTv\tNon-synonymousSNPs\tNon-synonymousTiTv\tNovel_non-synonymousSNPs\tNovel_non-synonymousTiTv\n";
	print OUT join "\t",
	  (
		"AllChr",
		$ref->{'snp'}->{'allchr'}->{'overall'}->{'number'},
		$ref->{'snp'}->{'allchr'}->{'overall'}->{'titv'},
		$ref->{'snp'}->{'allchr'}->{'novel'}->{'number'},
		$ref->{'snp'}->{'allchr'}->{'novel'}->{'titv'},
		$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'number'},
		$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'titv'},
		$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'number'},
		$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'titv'}
	  );
	print OUT "\n";
	foreach my $chr ( sort keys %{ $ref->{'snp'}->{'perchr'} } ) {
		print OUT join "\t",
		  (
			$chr,
			$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'number'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'titv'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'number'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'titv'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}->{'number'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}->{'titv'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
			  ->{'number'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
			  ->{'titv'}
		  );
		print OUT "\n";
	}
	close OUT;

	#3) Write out each sample's snp information
	open( OUT, ">data/sample_info.txt" );
	print OUT
">#Samples\tChromosome\tOverallSNPs\tOverallTiTv\tNovelSNPs\tNovelTiTv\tNon-synonymousSNPs\tNon-synonymousTiTv\tNovel_non-synonymousSNPs\tNovel_non-synonymousTiTv\n";
	foreach my $sample ( sort keys %{ $ref->{'samples'} } ) {
		print OUT join "\t",
		  (
			$sample,
			"AllChr",
			$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'number'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'titv'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'number'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'titv'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
			  ->{'number'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
			  ->{'titv'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel_non-synonymous'}
			  ->{'number'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel_non-synonymous'}
			  ->{'titv'}
		  );
		print OUT "\n";
	}

	foreach my $chr ( sort keys %{ $ref->{'snp'}->{'perchr'} } ) {
		foreach my $sample ( sort keys %{ $ref->{'samples'} } ) {
			print OUT join "\t",
			  (
				$sample,
				$chr,
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
				  ->{'number'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
				  ->{'titv'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
				  ->{'number'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
				  ->{'titv'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'non-synonymous'}->{'number'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'non-synonymous'}->{'titv'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'novel_non-synonymous'}->{'number'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'novel_non-synonymous'}->{'titv'}
			  );
			print OUT "\n";
		}
	}

}

sub writeHtml {
	my ($ref) = @_;
	open( OUT, ">vcf_report.html" ) or die $!;
	print OUT <<HERE;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Strict//EN">
<html>
<head><title> $inputVcfFile vcf Report</title>
HERE
	print OUT <DATA>;
	my $date = scalar(localtime);
	print OUT <<HERE;
</head>
<body>
	<div class="header">
		<div id="header_title">VCF Report</div>
		<div id="header_filename">
			$date<br />
			$inputVcfFile
		</div>
	</div>
	
	<div class="summary">
		<h2>Summary</h2>
		<ol>
			<li> <a href="#M00">Running command</a></li>
			<li> <a href="#M0">Basic Statistics</a></li>
			<li><a href="#M1">SNP Information</a></li>
			<li><a href="#M2">Samples' SNPs Information</a></li>
		</ol>
	</div>
	
	<div class="main">
		<div class="module"><h2 id="M00"> Running command</h2>
			<pre><code>$commandline</code></pre>
		</div>
	
		<div class="module"><h2 id="M0"> Basic Statistics</h2>
		   <br/><a href='data/basic_info.txt'>Download table </a>
		   <p><strong>Description:</strong> 'Overall' means total reads in the vcf file, while 'others' means records that represent other variants (e.g, indels) but not SNPs.</p>
			<table>
				<tr>
					<th>Types</th>
					<th>TotalVariants</th>
					<th>FilteredOut</th>
				</tr>
				<tr>
					<td>Overalll</td>
					<td>$ref->{'basic'}->{'overall'}->{'total'}</td>
					<td>$ref->{'basic'}->{'overall'}->{'filtered'}</td>
				</tr>
				<tr>
					<td>SNPs</td>
					<td>$ref->{'basic'}->{'snp'}->{'total'}</td>
					<td>$ref->{'basic'}->{'snp'}->{'filtered'}</td>
				</tr>
				<tr>
					<td>Others</td>
					<td>$ref->{'basic'}->{'others'}->{'total'}</td>
					<td>$ref->{'basic'}->{'others'}->{'filtered'}</td>
				</tr>
			</table>
	</div>
	
	<div class="module"><h2 id="M1">SNP Information</h2>
HERE
	my @header1 = (
		"Types",              "OverallSNPs",
		"OverallTiTv",        "NovelSNPs",
		"NovelTiTv",          "Non-synonymousSNPs",
		"Non-synonymousTiTv", "Novel_non-synonymousSNPs",
		"Novel_non-synonymousTiTv"
	);
	my @data1 = ();
	push @data1,
	  [
		"AllChr",
		$ref->{'snp'}->{'allchr'}->{'overall'}->{'number'},
		$ref->{'snp'}->{'allchr'}->{'overall'}->{'titv'},
		$ref->{'snp'}->{'allchr'}->{'novel'}->{'number'},
		$ref->{'snp'}->{'allchr'}->{'novel'}->{'titv'},
		$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'number'},
		$ref->{'snp'}->{'allchr'}->{'non-synonymous'}->{'titv'},
		$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'number'},
		$ref->{'snp'}->{'allchr'}->{'novel_non-synonymous'}->{'titv'}
	  ];

	foreach my $chr ( sort keys %{ $ref->{'snp'}->{'perchr'} } ) {
		push @data1,
		  [
			$chr,
			$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'number'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'overall'}->{'titv'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'number'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel'}->{'titv'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}->{'number'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'non-synonymous'}->{'titv'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
			  ->{'number'},
			$ref->{'snp'}->{'perchr'}->{$chr}->{'novel_non-synonymous'}
			  ->{'titv'}
		  ];
	}

	print OUT makeDownload("data/snp_info.txt");
	print OUT makeTable( \@header1, \@data1 );

	print OUT <<HERE;

</div>
<div class="module"><h2 id="M2">Samples' SNPs Information</h2>
HERE
	my @header2 = (
		"Samples",                  "Chromosome",
		"OverallSNPs",              "OverallTiTv",
		"NovelSNPs",                "NovelTiTv",
		"Non-synonymousSNPs",       "Non-synonymousTiTv",
		"Novel_non-synonymousSNPs", "Novel_non-synonymousTiTv"
	);
	my @data2 = ();

	foreach my $sample ( sort keys %{ $ref->{'samples'} } ) {
		push @data2,
		  [
			$sample,
			"AllChr",
			$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'number'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'overall'}->{'titv'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'number'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel'}->{'titv'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
			  ->{'number'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'non-synonymous'}
			  ->{'titv'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel_non-synonymous'}
			  ->{'number'},
			$ref->{'samples'}->{$sample}->{'allchr'}->{'novel_non-synonymous'}
			  ->{'titv'}
		  ];
	}

	foreach my $chr ( sort keys %{ $ref->{'snp'}->{'perchr'} } ) {
		foreach my $sample ( sort keys %{ $ref->{'samples'} } ) {
			push @data2,
			  [
				$sample,
				$chr,
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
				  ->{'number'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'overall'}
				  ->{'titv'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
				  ->{'number'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}->{'novel'}
				  ->{'titv'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'non-synonymous'}->{'number'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'non-synonymous'}->{'titv'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'novel_non-synonymous'}->{'number'},
				$ref->{'samples'}->{$sample}->{'perchr'}->{$chr}
				  ->{'novel_non-synonymous'}->{'titv'}
			  ];
		}
	}

	print OUT makeDownload("data/sample_info.txt");
	print OUT makeTable( \@header2, \@data2 );
	print OUT <<HERE;
</div>

</div>
</div>
HERE

	print OUT makeFooter();
	print OUT <<HERE;
</body>
</html>
HERE

}

sub makeFooter {
	return <<FOOT;
<div class="footer" style="display:inline-block; vertical-align:middle">
Developped by <a href="mailto:riverlee2008\@gmail.com">Jiang (River) Li</a>&nbsp;&nbsp;
<img style="vertical-align:middle"  width='35' height='35' src="data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEASABIAAD/4QCwRXhpZgAASUkqAAgAAAAFABoBBQABAAAASgAAABsBBQABAAAAUgAAACgBAwABAAAAAgAAADEBAgAMAAAAWgAAAGmHBAABAAAAZgAAAAAAAABIAAAAAQAAAEgAAAABAAAAR0lNUCAyLjYuMTIABQAAkAcABAAAADAyMjAAoAcABAAAADAxMDABoAMAAQAAAP//AAACoAQAAQAAAFAAAAADoAQAAQAAAFAAAAAAAAAA/9sAQwABAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB/9sAQwEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB/8AAEQgAUABQAwERAAIRAQMRAf/EABwAAAICAgMAAAAAAAAAAAAAAAkKBwgGCwMEBf/EAC8QAAICAwEAAgEDAwQABwAAAAQFAwYBAgcICRESExQhAAoVFhcxQRgiJTJRYYH/xAAeAQABBAMBAQEAAAAAAAAAAAAFBAYHCAACAwkBCv/EADURAAICAQMCBQMDAwMEAwAAAAECAxEEBRIhADEGEyJBUQcyYQgUcUKBkRUjJDNSsfByofH/2gAMAwEAAhEDEQA/AI7/ALhtk+7n8gbevArvVl1I5zR1XP6pQuIXdeDXK8DVBOY9Xtl4NTS0RrFtJJcL1TGpOC2UcITuiVJuewA3FVZUuTCxsZo0UkFnjkm5AshX2ihRtmUg8mrB/PTSzs7JhzWfciLFJHj2UZ7Mis92OVUFaJXkbwAea6LH8c/yue2uy8k4/ULjThbe6rHujy15wtPaz6axgN6v597hybpnTTLhaK6uh/wdA6ZU0qGjQXVwrYFosMHA+cCRRM157cbkQRJLviX71kpifuKkLannjnabPBBoUL6LDIbICiV+FKKVpSbdZDbXbKpqxQBAo3fBaymnxqPLvttLjeEWWfTG2n452+4JdtJPw2kznOc/W8WIPw+/z0zttHj8NsYG7D5nP3FgK9gTx+ePewOO/PReP7BRtaAvvRBHf+3xxyOqE+GfUFx9GMfZQF63RwkefPc3oHgNeiVKtluRef8APcVOWmkM5p55x2ribRq2IYMIv0JpdNxo/wBpD+31zJvMjgiQA0FSzZ9hx7kj3967fI64xslGO/6ro/8Ai/8Ax/8AvRBdDhNcR67kR67/AJaw6a7f++SXb9TGNddc512222xFJt+OuM5xprtJtjXTX8sJyGlsV+WI/ni+/uf4/j33tY+br47fzx7+w+f/ALPXlQ2urFrGDcKxozVaohoIzYitV5QC0lHuRo6HYFjlSDhTKJRSomkRM0Ei6UabUzA+8O+MdNmz2o9j+a/kfx8dh193iX44qqPzyP4+Pf36HN1/5evji4Y8Jrt69Vc83fjshVxCqmzNOiEhzsTzFZZRY1AUWDAoaglUT/mGU8n7RPDJgpiRANNFJv8AessAhee3H/v9urA+avdXkD1jqw286+i+RdbLXESDnp6bclTCxQbj5MjzPLXpphrDgWWIOeaEqVXpDMNDsZFJkPaHOM626t/pNpvr+eudc4zj7/5/6z/OP+v4+8fWcff8/WfvP1/WdZ1yZzjOuc4zjOPrP84z9/8AX9Z1nSVXkTohnrf5GPRHqjk9juFxR9XX9yc8pCvAFhqlhF5dznqPgtEiBXprILG7GR2FUCzsNMCKmEE3heh7/gk1LImmd2PFiae+MMlmX/iOjBB5lMHP2ndWxqO09z7cnqM9clzdSmlbCjjK/uoZL37WCNGp9QVNvmAuoZQWKEbWIZWVYX8/fJfWfEPyOWrjfLaZUT+X+o+k8IK6I8nELmjqGUXYPUPI7JZ64vSmRSHON6guoq8DMgJcmodej0n1kWwLBVqGdMeXHx3hLb0RzIWtTTTSbQFsk+lBfpqwBuJJ6PQHNJmLFfXKiKrcgeXjQ2d1WrBpH9I9Telht2N5hhrJ80bajfHjn0Fc0lAr3pJN67x43tfLWmScoFXRVndcVO4mhBLrUTZmqdVxcB50odmvMJHPxBK2hH2q8WhOwwwAv9vBoki1+RRsc82aFfyL6OCd0jIY8hCAbBFhQ1k89yAK7eogdgCraT8o/spEb2hX5NvWa5aO3+o+kejGJNJURpC+0Qeirir5FzmjJZXrOAmsAs5xFtl0epoA7pT1m7iyJ7WmaxDE6lFxFlxWdfujKKQSbO5tqHhqAsgMa7AGu56b8uoHHy0MzEJMZB5imthjjaWQ8DncisqiwC7AAkkdSt2X5PvSjbxTxHxB24e0L3xPVOZ2RZ3vgvb7Kc2K5XWLruIERVOm6O+o2DpD2e3XTlL2lv1V7e5cCkSqLdYC5VgcjfSLCTFdxMoIEXmEqxIKlgo27e5Ng0eKo7qJHWZuqkyYsSmRZcmXZArxhtrBJZAJbXbHSxPZIBtbCEWUMVxoxzSv7dm5f74WSicyuF8bdWYjQprbjoa7p72zdjNt5lNPvJjjpx93ufRHUNoWX+1xvWFjGeEWxiQwrkao0laKd/8AquYpEEO4OrJRBFlSpJsqR6VJoHsfUD049JkEyFxuKBmALqVYqt87SABZVgAKAJrhaPSHHf6yyunTTrcUBdBq9DuRPunSWOqHjabykyyk/wCnw1RyKdXBDvjbCQUzTOZ4oRtJzv20Mm8/CUCM43+5FIcsXGkb73Ulgu2QAAI1nkE/4PBXspBL0GHP2+pgPfihwOefwefbqvVX6jd+WdNhvVSIvFA/wrOBsnZBP2IdopLcXMhQLFHaxWY7/wDUCJjE21mw1yMX9YgPjI3GCG02eJkIBoszhAgNtZbbZXvQPJPsOtolMtkDYALJk9NC/fuR89ufaz1tpvgl9Kd19Y/GrwntnoS4VzoN5smttXjXdMaCQ8sdbrtnYIUpfRw1axWsT9NEyAaotoAQ2upW60CwT4gNflgh/GjZWZWoFWKmz7j46+LTfayn+5/x27/jowv6udcZx9Z/+8fz9fzj/jX7x/x/394/+c/f8/8AGpBABPvXHNiwTyK47f5465q4ZmSiCp7sKU/kH36SgrdT7HWPWvp3g3ljrhHNm/Dl/wAjyTnprVylSg8t5fZuoeVui1+lV9myX5S1yqVewXFoPTq+rihgq6nVdrPusjj2iwWVJp5IDICbhN0S1sXHC3fBsgWf4Y+m2pnT4kPmiFlVllQMVpNqpwWbkEm1BL1yRffd0BF5yYHhXyP8JrxDUzWr13iXig5k0uD5IYYhfr68O06LvajQZcQrTs3PFwsLCR5EmP1Fdf6gZxajshWBS14xDixGqLDZfJsCWUj8HuQKF3d9h1vFkpIsm0guszPtXaKDwQUaHI5FDcDdGvfrJuz9M63ce2VrFLsGekcDe+o+ndzptKa/4hzzi4+gVPf+u1KtG/5KAjRnZxemcvOrQNfxHYg1LRA4inHjXSCQuYtBjExiQCvR/mrN9rAoX8EH57avlrv8pjThwG78Wq2CSeCpNEAce/4lXz1Tlz/11nB/Pmt4B19NBVBrXkM7EQ9gJrWPa7XGgG+q5wqgOQNYtLKNKvmV2LQ5dAXFnfGgLZGsxEJxyeadzXwTFzx81fPwbviug+XtmzJ4WIJxYIpwpK8ecaBIIPDbDyQ18iiCR1ExXMz+IEeQ6K2p/TkE1D6f1wSpadHrxSp9eRqm1qnYBSAQ1cy8aGuK244ogP7U+B/rCligjNHJnjngUwMkjWeWceUARQJWzXFAdiSFI71yF6SZsDsdPzVsCAhyxZn22vl2XfcS53NZZqJJJBrphD+1J9V9Z9FeDe78e6jKrYU3y1aIKhy5lDW5k5MVYv8AX7ReHqNhZITihHU6xsya6jhDI17BOmZjwlMmI86wdY0dQbLkyxPMQYtiiRFoCRQi2xCjawVQigggAgpZIA6kqKOotkIpn4UfktQvubLHv7mxzdAINZ8bOPcfiBJyTkvn3k9Xzz3PDeic69G1MEesPbi3YB3uP17zvod8USunVouhHRQczqq1ZK+uW1RUlpOK3gMVpFqxq59Wv1OfSn6MahBjeI9Td8/LyDA+Hjxyz5+IQIWheWEFCsM6ygwujncqs3YRM8jeHPpzr+uYwnxYmG7YwZ9oVwxIIDMSdy7SpsAG1Hex1lnyw8zS1PyR73qdo8QcL808o4x2Xzkr8FdYoVRr9T6N1MB2+ErlirdjaJZWL/pv+SoC6zW282O0l66p7CKoAig1ZEQ/09/pZ9cvpz9VY8SLw74pwMzJyMTKzm0xpAmpYv7beXWeKRRKu1adiBIqR+tmClehviXwjrPhqVYsvGkp5Am9F9MhJAUAg8kluAauwACeSyD/AGnTXM3xB82Bktbmyzh9d7pBqucKpFQtR0/1YAdrWEk0ud9mSXOjDW2xs9M5h3Y2pstj/DdRNBBL4yRLd7aUgAiUSEtV+trtWANEN3rd/UB008jEmx0WfypI1ZVJV1dbDbiGAYBSBtIDWe9d+mY8S7412/jOv1j7x94/n/r7z95xnGcfz/Gf5x/Gc/ef5+lSIj1us9+5Io9xYDcHt8d79+uOT6MYzRkFuL7A9uT7fg1RHIocjrWe+90fpj49+m/KlVenuFLC4dPrNw6hQr8jhwGuttD9R+wuI2g86HfGgUitoTX6XZKIwTiyEnoC07jIZm4E4hU7mw5lKB6UOh2jdVBSQaAIrhLRRx/3USOWLqWm+dqMUUwk8h5ERwvpIEeO1bTRJDMu5iuwtuYAqCKoj6VvMvYHSj0TkvJU3ZvD785voKrKjWzv+f1LrlbZQazFaySm7qBbCmhPYzk7x4lFizLJkU3QTVLmyGQiNSpVY+AOR6WZgTZI/qbcaHFX2I6JafBDCsUTh93nStZZgxao1UEAEgVECAOQxIXksTJXNOk4S/EaJdl+ix8/4bcKp0NXEwhlEHFYsSIeYJ6swaRbasCogHFfmsIKOI0deQzPn3HEwK1z+57q8n7ZQoWhHXYgcWoHfjtx89+AaCc4yDOKF2EpmJLAWBG1gmmB5oiyebtzZ567fxo9Zap+jc1pIUKgex2jlAV8SOHkhy5fhtjk3R+SmsXkimXQjC0TTtq8mT/Hr5G4kc75mIWQGFuH/SnDnCY6wvW//edOwtWVVc2SBtBZQDfzQY2Okup4bJl5OZD3yI8bDyD98alN7oqKdrHhXJNMCCqt33dD74fyO98V9Jc+16HeUdlsP+put1Bgtr9zj6CuVG1ShwsdSMWlI1KUtNWmjGKU+NQVsbBAHHu40hZSmLF48VjxjaxDhncDeWtggHJNUBQo1ySSRZrotlHHmw/KCMIn2AkRCIjex7ItgAFgqqCVWgFoGunVv7WCqh8I+H3t3fGeDB5L/wBj69fphLSvxivGpucUeqU9ces0HBwY0QsZUhADaf8AIzTdnlqODmHSLeLVrZWakGK804YuqSgFuY2PICsvJ20oAFc1RpT6Xx4a06XUdbwtPJ/2cqaCCIIdktyybXpmNWgdyptRSJwfUW8DzJ25Jxaihct2VqgaqjkcbVWZDWi5YVexT8tuwTMVlTwL+QCsptrAkcEaQbbKogIGGf3cO5snh1+qz9Lfjjxr9SsnxjpORqGpY2dlCaSFckXHIW3b0aXftVbO4m7Ar0qSvXo1o3hibRNEwosZUSGLy1lOxWd5GW0UN/Tyt2RRAskgEdRn67LpPq6lY5PdxGr6lSWyqung5GmimN8PULOveJazolXGFwCo2DiXTDKaM4t4xCX7Ba2eHT9wAFNf6Uv0/a19Ntai8XZbx42rRYmVjJC6vMypmxCHIBdakLlLMZDEIWsKaBKnUPCmk6wI5dfhRApSQeYfMDspBNhCdqkRjcq7Q3IK0zXeT4c+8845p6CO8kw2iqoz7rzX/NVDmA5QwjBJ/tPAlgK/a14WbSJSNIifxawaTCQzMx1xEo2xcalnuJ6VeGJJmhdJYnQNNIzNIjCQvusvbAFi247ruzZPv1Wv6/YXh2HHxsbTFhx5caFVEWO0QjkjePYu5FIo2kUijgmgSCHotGQaZxHH97Zz9abY3z95z9Zxj7+s5/j/AI/jX+cf84x/+vxiUUbe9jk89/fivb46qzGpXDVZTQ28ixdDnuSDZocWKPHbpOT+7G8N267eaHXsnlI9a2m56loiP0PlyG0b3P8A26r18X786/26hhBMERrltpv9ofdP3mnCXnp0NdNk3iMTGSMthlsnpUmj/N7h2+bNix8Ef4yaBJsgOd1rNQq6I8soQ/H2UxJ9jYB7nqsnx7fBz6gNbeFx/QKXntx8wMfBXXm9sxDzmmh3PnnV+wLExFd5FcS31ybwPXoKCwrSknS1VQ0EHiTW6ovU5Mcay1m7rlGwWHoA20R7njkHuKsfFfgkdJZtPRpkmBbzN7uGDWtAsaFkkEtQuwVskV2K+/vbxL75+MPwRxHkPoGv83q3NfTfUH75qrVL92/U6tZeabJ2SyjdFtq0wqnbJnQr1leEyRVlttOSuYz7tCoYoo4Fv7wiIJYs0Cew73dfaSRZJ4Fkj8n6MGM5BzrYttCUSCQCtDaAN3HbcT+Se56YJtXwpduF+PD47vanhlRVa16t4R5xr/Yu3czu1HO0unoEl9w+ltyOWwJVySY+ZhqNHZKkbQWZ6BXdN3mGppId5OKN2R/vXsdyAoHpLHspsWARQJvtwFHzYybCQwPvP3u8l9yCXsAAE3QtRYPB+2uD7P8AbcB8++R+v+z/APxQ+fuOQG8g7Bw1pVq5X+fboSUVgYY6++sLB0zsjCxXt2wbMK0vTugLg6NX7KUo9dmVDgxnQ7Y2TIeKN0aBBPfsBZo2RR/uaN11qdLSSEICbFGy3IIYEkkGvgADsVoE1zBveL92jy9p4b/t7vOB7oa/2T0VbupeoBeMOHh72ucIuHbH3TuY8cRXKcoONhOXyKOKydJNPD3jgjEQqWhsME1m2bj8xJcjAXHy0QM422otbG8LXNWLuwQpIurAHRvSM/I0fOizAylkmjkxmPO3aw2gkDgkKwJBB2MVvmzYD3BPzbzx7WI850I0i0who4L5YE8ZUz5/x5c6rh7kKqXGYeFiUTao10w9i1k2yZBvWnS6Sb9AyFfMZBnjCB9N1mDBxd8olVn2A7rQsV7sBZAAsgiw328sBev6cfVKXVtJxMXUGRWjVY8iV12JNIoHriVnYpuYlE3SMQVIJIIDcXIyef7XxDJ2lBczuZ5N0GzV6g8EWW63GYwvmTiMzJ8xFokBceCcMh0LYe3TafiLu4Rx5YB6HdJWDFeL/ieTkPxueiNxrkrHTBRZPsTRFkdH/qNqGoPoTzaLMMSQoW3zkyEodu4wiPcvmKCDH5gePcP91XjtWUn9NLO2fGT8m1y6Jx+2t9LBxDtK/rPCr6/zO9aNecNiMWPl1lswjqQtg1CcU5lBU7qC3yRtlhraam8yKf8AcM0qRw7YozvjkZkUs8ZtXYoCQwoVxRo1tIK81fVANVy9Sy8udtTklknEhVzK+9iULItsQCwC8KeeCT3LXtZ/jR93UT5GfHPI/T1FlhEltqWJffqvEZDKRQ+n1nAyq+UgyKLecyAZW803NTTlZkma1Jmiebybjt4JZFaEGuaUDjmhYoHk13s3Y7H46DG5iYewo2TyD3PexV3/ADXH56u/Z6tX7ahOrtsSrLKgbafpGp7CpXuVpOk22ukY56hqOUAaNiXEcmw5Q0sOc6fe0ec40+uAUHJ+0MpVj2Ff/KiL7A+wPt79dt5DS0eDY4+CfY8VVcHj3Hv0Av4f/WLjqnVvk2iylbFUYX1X0frfPsOX1e/zyquvHpdHiom9cjsTQxUqUkUId4ExRLZKcCxstq1ZPpXxOocq7MhIMa9yAtmgbbb8D3JJoj3BHPQfSnBGZVbRkOR34FsSLNmhfyfxYA6gj0N7j+L/AOY1X488a3ah9V7ZV/Rdpot5dueXnwqFfm7oapJMauonUbpEJpMqeEMGT6i3NTVYp9lisqORucrms9ZXH5+2Yx7q3BVayPUbUAkGj7D5Fe1j3W406DzFLEEexBFFmJFDjjiwRwwNrYN9HFd+rfKnE33AuFkdJqiRh1qOOmchR1uCRoikEruhlbBHiPqwp6SurYnKLSlrjmRwC3/VuoNU0nhbThC/1wXEylSR8dhuQ26naABV2Uava6aqPz8/ElQsAW2lm4BIogmgABx8g1xzz8Dk5p588reG6x6L6fzSn1Xj1e6Jc7x6U7vasb7jrDrJMmw4t1saSGbbDqEatQuNOHWhzhpEX7txKvBHIYG7lpRPJOxhMbGZRt3U33Dg0STXNc9gOR89K1s/b2NDi6783wPz89K7/wBuqlj95e6/kZ+XDoKUyGzuumFcx4wE7G0mMpYV8XiumcIksEuBMFpuUBcmoMJosQWf8eRZy9psZblY3xpisQwM2cpksrFDuYDkt5YUtZDWQVFGxVUxvrJUjKGNOG3g9u3BLGl7XbE/gm/fqWPjTO557F+aX216mm6VrZF1KbXlbyNe0WQhFvaovmVcYVkI5xn8eo9ZqwqBnMOqaIdibEFaAHLWPSfSHYWN8HTpc/xPNk5+WuTjYcUsUcIouWZw8bDeCfQpCqBRKoPk3LWt5OVpXhLRcYqFE0kYkItbv1tvKkSMXK1e4FWIIPA6wf5F6Id506Z2y+KOWPM8wdXh2sWaJntfXat0Davp5J56/I+ONHT11M2OsIse8aHZYqa4mmg0a7j/ALQcZ4r8RaP4blXMmxnXHabyQT9th9gWy1d/i+/+XboHjmfJwMXAlkaWbHxRHLuO47ttcgk88DcxssSWJ5siv/uWKNy30Pzfxj7982oILfTc8tk5/wChWKoUWZzWg5H9aZpyL5CDgAyFuvsTp7S7Y2jhlWqLS1SCtJwTrekkPfnhbVtP8QadFkYDERs7WA1ksB6h2NkDkkAglrBI6gnVIZhqWZLMrLvyJJKO4Ah2J3CyeCLAA7VQ9+pE/syvWTJF2n1B4nYHQSVO81cPvVBEPm0gm0stGPUUq56qhoVsmxRDumtq6c0xMzHgD1qA0kcMxM22+rim9P8AsgjfZVq72bPNc+oLYPNHtRJse9NuoAX2A4Ht/FX/AO1362FrkqQNUedELOfKILKbCCL/ACSfMHFuVACPj62+5zJodBINc6Z/KSXTX6znP9dhRyFHpHoZb44ugS3vQDE0KPevcdcpD5LS7gaAPbltwHIAF9wRVe5/t1ri/C/ySxeF/V/QL50flnQmNU77wv1Fve6NXEyczo1TtFZ6et6BRoG81sa1QSuK1Ld/Y6jdTJDNYR/8hC2JA1hVRGamZoWnyYYkKhiA3rO0HbQsmu5FV7knuBz03NMyVXF1CU3ayN2U3ubzKFC/+0Xts9wDu6qX5mq1U5f6D5H0ikprqj7q0ZddvPnexTuUEdertU83WbNHhou9FcV/Wd5vcL1Sejm2xs4LyaK6/wASKKmyYeTbdi8OIip6nibdI0K0QezBCpHpJVnstwALAYA89NSTWpIMmSIw5YWFMd5ZBG7RM+UrvGXkoKDs2pECaCLtWzuuGeLR1uzdW43eLhvdA+w986Z0MnhzmuXnKqjcLaqEarsd/MW1c+Yo/IFu7F1B9SR141ir8YANciZZzbbMY3kjFa7hDExIc7GkmbOmyGxxAguMJFs9RPIN7yOSD2sEdGcPJOTkNFHkIghVJnkmYIree0gijQ7yQVig3WqglnYFVAUl+v02Oz99/EFBQnds2470H2H5654snm/I58zTuHeajYLnpqCK0GMbQRhrXy0+LRnKv3hOwI2lYg5l1Irl9dvrron0O8ORa7qEbZmoPOI4tOx3WOeWUBWKsyxyOsLorqWEbMLJUBiNsu+D/DeoeJ5kxcGMHam+TIdT5SoAVJ3EUHLMGQMSLo0QCOqZ/CVSKV8dcPQPjZ6F0HUTqXUb/Zuy+eb2eJFWA+yVUmoVNA7VoIdmTCELqHPSEohtwoWzslg0RsVd6qOrCvkto6+H+iX180b6/wDhka1p+NLh6xpWQmPqWnzyBp8Od1eSNhMEhklx5UG7HleGKyHx9rvCzyFfFvg7K8J50T5JieDItleJwySbAFaMrQANghlAYCxYAIPUm/BB4L6H5YrHY+p9RPrttsnbGboyq3EBic2Z2RGt6BdJS7a3LZbbSwmdJPL1umkGfyawrmIItk/SaA5FHljQ9PnxMzM1NxD5k7h4lLnYIAdjlO2w8ABALItbomufifX1zodFwyzDHhpWDnkuNtbhRut1bjXfsABc9fKLx+S2cxnQAVLWynm3OKYTG0ed9hak5CIPcQfufzjkBIweBrHqbDNoSDtLrKJLrP8Ap/hXb9Xni3H8K/TvCnkijhyn1WGXzlUBPLbL3Lbdg1bQV7fnuoP/AEm0/wD1nxTkaeULNJEpQUSlBGJAN0avsAR7g2a6Ty+S3pvQPH/kTzJQxanX2dQTXv1dSCafadpz0VgqndeO1usW+q2UDQcf/MaIHybe8DshjpzKvfUVYZRQ4IVjsRm/+kj6t4fjjG1TTI8hsjO0tmy8gIlxpFNNBEjI6gKJJwyBQ6qrIWKhpNqMc+r3hNtCy49RhhRcZ1GIY4mBlbIiRt7eWRuMZvdYNs3APHGaf2iPmO3svZ1g9ghjmqeTV7kfTuUTMLDCNqTZej2gimH5QI94TdCGY6xCLvYzXOVkgsc2CFG+cMpdP29ldS+pPhbTPHWm+Bs/UkXXtYiTKwsNCr3Ef3KbnfzD5crukm2LYWZVSRQEKs8TR6Fqb6fJqn7dxiRnaxKsHs7ey7eV5FtYAJogmgNkaRrjTTbG+2u2MTC64/L7xjP1PHtjGcfW2c/+bGP4x/zn6+8fjjOcyGnqnDC/+kzEAX7i+7D2HA9vz0DHreZmIA7ih/J4snmuB8c+5HWoV9vNAqr6Y5lYrSxZp6ywvXrivHs04jUkgkQ98apMXR6o3aCYgKxiWkBWzHweKu2XtyJLAI0RGMVJrolgLZOIBvDMFbdyDQaMEcgH+ocAKaNEAX01MWWE4GqpG6blmayzBa4dgbJA4Ckr7kjg3x0Q24en/Ovmfj3j/WzdI9CRnNC+59fR1DmNToFmC0rtg9EdUgV0li7s1/rhgYLQpCOzL2SClr52M+50oJRP7cKfs+O6Y97jZybG3uak9u1g7QKNCyB7WA4y/wB5mTwpjxhVxcQNPJaqLTbYcIQTQ3kHkKGewA1Zf1T0lwbg/K+Mtur3nr1fuPWO2dH6olR8MoXNrkq1wgxQUdiq+Xb64UdsuhtzI0qYLYUcdedHLLvMNpmCWD+kWphMTDd8jLMSRPmSoCtqGij3ne1+kWRyBfB+OlWm47ZeccODTUyGjXHZ8neA588SMiohUtIkZjO5xwm6Lds3DcXf4fbWqtdI6OKNqdphDfrbAqwyn3KNASWS0F2LYLSGUk4JTtEWZoQYMikFXTsI9SZNCiBxC4fz4/rUm8Q6947jzsnXZv8AS1ymyIsFpG2yI6uikLvAG4AkBVPFfN9X3+ksGMnhhsNcUY2QGG+TaNzcqff8kXRPubAIAsj8mPkFh6282vEdIcH03uXODgeqefL3XjiEj2pdZpWuzCubrnwcsbBVh5rqWgIOEl1JE3YAMRcwkKh87Rf+mn6z5X0X+q+DreVOZPC2quukeJsMljF/puXtWbIWIMitJjuI8mMFikkkSI6MhZS7vHXhOHxVos2lRkDNx4v3sMpBLCSMFVj9VqFe6KhQaLEHtRQvhdsvsOz/AB6cNK9yU2Gk9tRDva/rBpChGPsHP0LHUOi2t+JXjzEw79yiwPow/ZFR4YzAZayghmHTxafoZ0vN0zXNO0jP0bNizdNz4o8nEy4S3ly47lWB9YVwG3AqHRHBADLzzR7UsWVc6SPJUwvizmMxtVq0ZKkN7EcLfYXdGjXU2fIDD+w4qPeogyioahYVGllgExDGTLULNg6stJpJZp4dQYlR7sF1qTPtrCLuBBNLJrF+cuK1/rH+nOf9RPprk6ZpwJysbLwnjIogASG9wPDAsAKFA3RsEjp8/TrxOfC/iSHUAoO5fLBFhh9qgAi/uFjkdzXF9LpdBoFT9lJ1vG+ncwB6PUGt0VtmYoqaJ+MsFLxJGXY2DaUmLapy6tBB36dslI2w6Xv267dtHGxjCI8vfC2N4i/T3ga/rXh3xgMPxHJg/tnwoplXIifE3MWmjMjPPHJYQRyoVUooO0op6sjmaunizVsKLP0hJdKcKwlkQNGrULUAqRY57HcBYFA9FZ8xc+oflCt0Op8jqCup02iC7hgpFkGgQeRT8klHRTDQx65ycToXOdIWR+vLObnMhMuN/wAttoT8HfXXxbjeO8f6t6vqOo614jM5lkGfLPIJCVYCOASu644VSQoRQqglVVV4BjxJ4M0nL0nMw8WZMGBkDBQn2BWjPB9KgFlUVYYg2B7qcJ4vHaKD1peZMDsB91xGYpMRSfoMtNgSMxyfpyfpyfoEyZ02/DbOkn1tj+cYz/X6aIL81ZDRPJoixwQ34uio+L9+OqJTKY0YLxvAD0eSpKq1cWDV8m647djrD+r+LlnbKTz4Gnc6qdojrHoy+VfnsbTt9LqjKXk7sjn1I/1QYkS3BG5cNArDzVuHBJKlHK3YQHxaV2fBa7chyPkmQROSNyspobrI4J7EWfSAt3VnkA10xo8BoGykjZ9uQx3emOQkgOFCmVHZAwkbcVpm9PNKldj0v41fN6jzmt2+qU516a85dK6PwHolmtvUFKKl1GwCXclrmZtSiB5E9lMsZFqDunJlwIClnYtb3WlkdSszSWCvg6HJLKoLGlfzALJUnduBrt7Aj27AgjjojHhyL5qCNB5saQyP5dSbY12GuAQwshqstXwxBxG4+ObjJwjnvO+xcsrSz1d486hLSLoivt9rtLApQZtZr9gV2a/hnxQj2TnVjpSWv9DTlSELRX0kziT/ACLGwTuEMwnXopdS0mSBSpcLkBWIskyRqGsm1o0pugRbAVbHpVpwfSNVOTG8qpOqEoTtRRuJXao5RkDMu3itxHZgSxB8f/m70z5b506t1x8odkXI+msSrQLokBptzdLK9NJEWpPZ02mWtha0WzEQoecVTKjLfj4yQOzSKyQTBhvGv9VX6b/qf4r1Y5vhnHiyIsZd0QjheSOSNSHNqzLIXCvt9SBS6ttsKCbi/TXx54cx8RMTUcgQu4Q27iMlmNglmGxjdFlRt1e9EXe/T0LzE4JkMRYw1zRdNMuap3/7hE9Vsof19TAWCRhDAxCMFm0jHlgk0xPpvPpLjGv6f9eaPiXwP9QNKy49Gy/CerR6lNOkcjeVIIyd4QMw2WQSt9+3qDdup501sHUMyWfB1TF/Zyr62Z0Ztn27UdnAAC2aPfsST3E8+7n/AHB0/Xugi/FpybS6eXGsMD1Ew7NRObJKektBhjHSzj8tt3T7LQpHq9kXtA8YrkpNjRrnZzacP/GzlMQdP0F/o90nxHp30V0LB8UR5SahjOwiOQ7mVcRkjdIwpSkCu0y0spBTafKRgzNTf6w6fg4vi3J/YZMEkUu1pVx2j2iQ71shXZgHRQd21eABuYHoaHrbxp/dk+nmCe998R361KucNQ7jV+cc7615vAqwT1NJiZUxW8mo93GW3FutnzruHq8WWZnLjaXTTeb9ciKWzmZiwaji5OLkRJNH5ZLIQTuCneQzqGcEqrLuQBhu9Is8x9ivDj5mLe5VVlZrpTuA++32cA8bWbbQu+K6Pj5V7rUKDwqroujtJU91qlfDFvcVlG0W2iS0rNMDvcvVG0cJ69+HL+7DZClDQMBJh9lRwQs8MMOPzdfWfwp48H1f8UE6XqGcmXqmXHjY+LHIEjgad/LG6mtEU0FsXtH56vVoOboreHMCTzY44oFJ9LhJL2V2AFgmrJqh2Fmui5efeTWvpsMd1vCBzRqRtIEVXVTbWIW324GcWM3Qp4pzrtmoVuaEkfXVM1gHs7aPBOWEFcH0Dwfc39Ov6Js7WMDH1z6jwJ/pEjHJg0uBf2ucm6LcUknLEBfWG2CMM43Na0paFfHv1UXDGRg6M8UiG0lyZEDgMCwPkhwod12EGQHZG1Iyu+5VKSVHj9DfG8mItc6/WZM412zr97a4+8fljOn395x9Z2x9Yz/Of5zjOPYhaUihwBwP7HqtbguOeeQT7dh2FdLKe1P7a7z16v8Ac9N9q0nrhHnRioe89t1p5pQeNUk6sXW7UOzauzLXI1gtNTkVPLUCAsVWIkVMwyVsogdSbFHSMf1ury8BhVgnmwTXJNWNwJJFndf9PYjrikITcTXPFUSSPVSkbtrVuJX02LPWYenP7fii+mflErfyZN/Rrqutqve+O3jPGI+R1N/WmE3HAkIUC4q1sL2ER/65uhF2JO/0tmZZNPvjEJG/1Pty849tx4/n2/t367Kq9woF9+B7H/B/H9us5+Sn4C/P3yJ+luH+r5r+94j13lxta3vTBFRlt3A7JXqU+VvKontSx0/SBAtEUkDhdHYYQ2xbBS2hCZAGDpgtJc887Si+okAdzfF0DftR/tzVWesJi3KhAvvXF2fjg+yntVAXfuD7ajw51/LaGHONvyzrj9OPONY99s7Yj11zjONNY8b401+s/ec4ztjWP8ttMcdxQkywobBG4hST7V87TYHfsee5PWjSGLmORgpIG2yADd2FAPI9h2oG+VA6h/o/nnh3XJwzOncpoF8MWxGQhF2msKXRQ8DD9HQyLUg0WWbaEjEEWN45Nt9Px0xjGun3n8mZq/gzwZruTHkatounfvUYskrwRs4JKm9xSySVF03pPIPRjE1rVsWIw4Go5EJeg5SZ1QgWQCqEgfcwFAWSAeKqVly8FYOMEAGKGMJBrCNALBFBEMPFHpDFBBFFprrDBFFrFHDBFjWOKKPXSPTXXTXX+nbhQYmm4sWHjBVhRkSONAQqj2r245scVYP9XQ4ySNIf3OS8+QdzndZsbixYEgA+om75JaybLHrv76xyZ/CTTTfG2P50300zjbGucZzjOPr8c4+86/xnGcZz9Zzj+fvO4P7dMmZgAosk+4A4H82rURftx351VWkbeLIQFQeTwSaBq/z79911VGBDfLfnJhetulHcU5aVet20b7e1T0evyut3cUedY3Ehu4O0kjLX6xJgyT8p8zYwRtJ+511kwzI/p74Nnzv9Wl0bDm1ScjIOS2PHJId53AncoFg8G7N+xBsF08SayMY6dFmZQxltNnnuqVypUAcgfAAFAD+Op1xBBFrtiKLXTXTX60xj+Pr8tvzznXH5Zx+W2347Z2/nbO2M/ecZz/Lyhx8ZEGJjqccIKGwVQBoGrr888dvngRNKu4SZjFge90QBVUKHPPFkX29z1//Z" alt="message" title="riverlee"/>
</div>
FOOT
}

sub makeDownload {
	my ($in) = @_;
	return "<br/><a href='$in'>Download table</a>\n";
}

sub makeTable {
	my ( $headerRef, $dataRef ) = @_;
	my $table = "<table>\n\t<tr>\n";
	foreach my $s ( @{$headerRef} ) {
		$table .= "\t\t<th>" . $s . "</th>\n";
	}

	$table .= "\t</tr>\n";

	foreach my $ref ( @{$dataRef} ) {
		$table .= "\t<tr>\n";
		foreach my $s ( @{$ref} ) {
			$table .= "\t\t<td>" . $s . "</td>\n";
		}
		$table .= "\t</tr>\n";
	}
	$table .= "</table>\n";

	return $table;

}

sub isSNP {
	my ( $ref, $alt ) = @_;
	if (   length($ref) == 1
		&& length($alt) == 1
		&& $ref =~ /[ATCGatcg]/i
		&& $alt =~ /[ATCGatcg]/i )
	{
		return 1;
	}
	else {
		return 0;
	}
}

#return 1 if it is a transition
sub isTi {
	my ( $a, $b ) = @_;
	$a = uc($a);
	$b = uc($b);
	( $a, $b ) = sort ( $a, $b );
	if (   ( $a eq "A" && $b eq "G" )
		|| ( $a eq "C" && $b eq "T" ) )
	{
		return 1;
	}
	else {
		return 0;
	}
}

sub info {
	my ($msg) = @_;
	print "[", scalar(localtime), "] $msg \n";
}

__DATA__
<style type="text/css">
 pre {	
   margin-top: 0;
   max-width: 95%;
   border: 1px solid #ccc;
}

pre code {
   display: block; padding: 0.5em;
}
 @media screen {
  div.summary {
    width: 18em;
    position:fixed;
    top: 3em;
    margin:1em 0 0 1em;
  }
  
  div.main {
    display:block;
    position:absolute;
    overflow:auto;
    height:auto;
    width:auto;
    top:4.5em;
    bottom:2.3em;
    left:18em;
    right:0;
    border-left: 1px solid #CCC;
    padding:0 0 0 1em;
    background-color: white;
    z-index:1;
  }
  
  div.header {
    background-color: #EEE;
    border:0;
    margin:0;
    padding: 0.5em;
    font-size: 200%;
    font-weight: bold;
    position:fixed;
    width:100%;
    top:0;
    left:0;
    z-index:2;
  }

  div.footer {
    background-color: #EEE;
    border:0;
    margin:0;
	padding:0.5em;
    height: 2em;
	overflow:hidden;
    font-size: 100%;
    font-weight: bold;
    position:fixed;
    bottom:0;
    width:100%;
    z-index:2;
    text-align: center;
  }
  
  img.indented {
    margin-left: 3em;
  }
 }
 
 @media print {
	img {
		max-width:100% !important;
		page-break-inside: avoid;
	}
	h2, h3 {
		page-break-after: avoid;
	}
	div.header {
      background-color: #FFF;
    }
	
 }
 
 body {    
  font-family: sans-serif;   
  color: #000;   
  background-color: #FFF;
  border: 0;
  margin: 0;
  padding: 0;
  }
  
  div.header {
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 200%;
  font-weight: bold;
  width:100%;
  }    
  
  #header_title {
  display:inline-block;
  float:left;
  clear:left;
  }
  #header_filename {
  display:inline-block;
  float:right;
  clear:right;
  font-size: 50%;
  margin-right:2em;
  text-align: right;
  }

  div.header h3 {
  font-size: 50%;
  margin-bottom: 0;
  }
  
  div.summary ul {
  padding-left:0;
  list-style-type:none;
  }FilteredOut
  
  div.summary ul li img {
  margin-bottom:-0.5em;
  margin-top:0.5em;
  }
	  
  div.main {
  background-color: white;
  }
      
  div.module {
  padding-bottom:1.5em;
  padding-top:1.5em;
  }
	  
  div.footer {
  background-color: #EEE;
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 100%;
  font-weight: bold;
  width:100%;
  text-align: center;
  }


  a {
  color: #000080;
  }

  a:hover {
  color: #800000;
  }
      
  h2 {
  color: #800000;
  padding-bottom: 0;
  margin-bottom: 0;
  clear:left;
  }

  table { 
  margin-left: 1em;FilteredOut
  text-align: center;
  }
  
  th { 
  text-align: center;
  background-color: #000080;
  color: #FFF;
  padding: 0.4em;
  }      
  
  td { 
  font-family: monospace; 
  text-align: left;
  background-color: #EEE;
  color: #000;
  padding: 0.4em;
  }

  img {
  padding-top: 0;
  margin-top: 0;
  border-top: 0;
  }

  
  p {
  padding-top: 0;
  margin-top: 0;
  }
  
</style>
