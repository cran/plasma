#! perl -w

### The goal of this script is to (a) take an input list of gene
### symbols, (b) convert them to Entrez Ids, (c), pass them to
### ToppGene to run a gene enrichment analysis, and (d) parse the
### results and write them into a tab-separated-values (TSV) file.

use strict;
use warnings;
use LWP::UserAgent;
use JSON;

my $json = JSON->new->allow_nonref;

my $ua = new LWP::UserAgent;
$ua->agent("tricolor/1.0 ");

my $httpHeads = HTTP::Headers->new;
$httpHeads->header('Content-Type' => 'text/json');

my $verbose = 1;
my $logfile = "toppgene2.log";
open(LOG, ">>$logfile") or die "Unable to write to 'logfile'.\n";
print LOG "START\n";

###################################################
### This informationn needs to be changed every time we run the script.
###
### outfile should tell us where we want the output, probably using the
### name of a microRNA.
### my $outfile = "myresults.tsv";
###
### symbols is the list of genes from one analysis (of a mir)
### my @symbols = ("TP53", "BRCA1", "BRCA2", "ERBB2", "ESR1");

### In the code below, also see 'faves' and 'catlist', which we might
### want to turn into user options some day.

my $infile = shift or die "You must supply an input file name!\n";
print STDERR "Opening '$infile'.\n";
my $outfile = shift;
print LOG "Using '$outfile'.\n\n";
my @parts = split /\//, $infile;
my $fname = $parts[$#parts];
$fname =~ s/txt/tsv/;
$outfile =  "TGOUT/$fname" unless $outfile;
print LOG "Still using '$outfile'.\n\n";

open(SRC, "<$infile") or die "Unable to read '$infile': $!";
my @symbols = ();
while (my $line = <SRC>) {
    chomp $line;
    push(@symbols, $line);
}
close(SRC);

print STDERR join(" ", @symbols), "\n";

###################################################
### First, we have to translate gene symbols to Entrez IDs
my $lookup = 'https://toppgene.cchmc.org/API/lookup';
my %hash = (Symbols => \@symbols);
my $firstparams = $json->encode(\%hash);
#print STDERR $firstparams, "\n\n" if $verbose;

### Set up the request.
my $firstRequest = HTTP::Request->new('POST', $lookup, $httpHeads,  $firstparams);
print STDERR "\n",  $firstRequest->as_string, "\n\n";
### Send it.
my $firstResponse = $ua->request($firstRequest);
### Test to make sure it worked as expected.
unless ($firstResponse->is_success) {
    die("Failed at the first step.\n  Status: ", $firstResponse->status_line, "\n");
} else {
    print LOG "Successfully submitted gene list.\n\n";
    print STDERR "Successfully submitted gene list.\n\n";
}
### At this point, we were successful and have received the HTML page
### with the JSON formatted response.
my $firstContent = $firstResponse->decoded_content;
print LOG "Still using '$outfile'.\n\n";

my @hfn = $firstResponse->header_field_names;
foreach my $fn (@hfn) {
    print LOG $fn, ": ", $firstResponse->header($fn), "\n";
    print STDERR $fn, ": ", $firstResponse->header($fn), "\n" if ($verbose);
}
print LOG "\nCompleted first headers\n\n";
print LOG $firstContent, "\n" if $verbose;

if (1) { # save intermediate results to use for second step.
    open(OUT, ">out.json");
    print OUT $firstContent;
    close(OUT);
}

###################################################
### extract the entrez ids
my $decode = $json->decode($firstContent);
my $garray = $decode->{Genes};
my @entrez = map {$_->{Entrez}} @$garray;
##### WARNING WARNING WARNING WARNING WARNING ########

### Do not *touch* the @entrez array until you are ready to
### json-encode it.  Really. Not even to print it to a log file, or
### for debugging, or to provide progress updates to the user. If you
### do any of those things, or ever look cross-wise at it in a way
### that might cause it to think it contains string, then the
### god-awful stupid perl JSON module will encode the integer Entrez
### IDs as though they were strings, and then the ToppGene web site
### will toss back a "500 Internal Server Error".
###################################################
my @official = map {$_->{OfficialSymbol}} @$garray;
my @entered = map {$_->{Submitted}} @$garray;
my %mixhash = ();
my @remove = ();
foreach my $i (0..$#official) {
    my $off = $official[$i];
    my $ent = $entered[$i];
    if($off ne $ent) {
	$mixhash{$ent} = $off;
	unshift(@remove, $i);
    }
}
foreach my $key (keys %mixhash) {
    print STDERR "$key => $mixhash{$key}\n";
}
print STDERR "Removing ", join(" ", @remove), "\n";
print LOG "Removing ", join(" ", @remove), "\n";
print STDERR "Started with $#entrez.\n";
foreach my $i (@remove) {
    splice(@entrez, $i, 1);
    print STDERR "After $i, now have $#entrez.\n";
}
print STDERR "Ended with $#entrez.\n";

print STDERR "Got the Entrez IDs\n" if $verbose;
print LOG "Got the Entrez IDs\n";

### get ready for the enrichment call to the API.
my $enrich = 'https://toppgene.cchmc.org/API/enrich';

### Put a 1 instead of a 0 for evey category that you want to see in
### the output.
my %faves = (GeneOntologyMolecularFunction => 1,
	     GeneOntologyBiologicalProcess => 1,
	     GeneOntologyCellularComponent => 1,
	     HumanPheno => 1,
	     MousePheno => 0,
	     Domain => 0,
	     Pathway => 1,
	     Pubmed => 0,
	     Interaction => 0,
	     Cytoband => 0,
	     TFBS => 0,
	     GeneFamily => 0,
	     Coexpression => 0,
	     CoexpressionAtlas => 0,
	     ToppCell => 0,
	     Computational => 0,
	     MicroRNA => 1, # may be a "positive control"
	     Drug => 1,
	     Disease => 1,
    );
### Create the lst of categories of interest. The parameter of most
### importance here is 'MaxResults', which determines how many
### enrichment results we are going to return from each category
### tested. I suspect that 10 is too small.
my @catlist = ();
foreach my $key (keys %faves) {
    if ($faves{$key}) {
	push(@catlist, {
	    "Type" => $key,
		"PValue" => 0.05,
		"MinGenes" => 1,
		"MaxGenes" => 2000,
		"MaxResults" => 40,
		"Correction" => "FDR"});
    }
}

my %parms = (Genes => \@entrez, Categories => \@catlist);
my $secondparams = $json->encode(\%parms);
### We should now be safe to touch the Entrez stuff.
print STDERR $secondparams, "\n";

### Set up the request.
my $secondRequest = HTTP::Request->new('POST', $enrich, $httpHeads,  $secondparams);
print STDERR $secondRequest->as_string, "\n" if $verbose;
print LOG $secondRequest->as_string, "\n\n" if $verbose;
### Send it.
my $secondResponse = $ua->request($secondRequest);
### Test to make sure it worked as expected.
unless ($secondResponse->is_success) {
    die("Failed at the second step.\n  Status: ", $secondResponse->status_line, "\n");
} else {
    print LOG "Successfully submitted gene list.\n\n";
    print STDERR "Successfully submitted gene list.\n\n";
}
### At this point, we were successful and have received the HTML page
### with the JSON formatted response.
my $secondContent = $secondResponse->decoded_content;

@hfn = $secondResponse->header_field_names;
foreach my $fn (@hfn) {
    print LOG $fn, ": ", $secondResponse->header($fn), "\n";
    print STDERR $fn, ": ", $secondResponse->header($fn), "\n" if ($verbose);
}
print LOG "\nCompleted second headers\n\n";
print LOG $secondContent, "\n" if $verbose;

if (1) { # save intermediate results to use for debugging third step.
    open(OUT, ">final.json");
    print OUT $secondContent;
    close(OUT);
}

###################################################
### Parse the JSON results and write them to a file

my $decoded = $json->decode($secondContent); # hash reference
### should only have one entry.
my $anno = $decoded->{Annotations};
### Should contain a single array reference
my @annos = (@$anno);
### Each one of the array entries is a hash-ref;

### These are the terms we are supposed to find in each entry/
my @terms = qw(Category ID Name PValue QValueFDRBY QValueFDRBH QValueBonferroni TotalGenes GenesInTerm GenesInQuery GenesInTermInQuery Source URL);

print STDERR "Terms: ", join(" ", @terms), "\n" if $verbose;
### Well, there is actually one more term, but it contains a list of
### genes instead of a scalar value.
my @colheads = @terms;
push(@colheads, "Genes");

print LOG "\nTrying '$outfile'\n\n";
open(TSV, ">$outfile") or die "Unable to create '$outfile'.\n";
print TSV join("\t", @colheads), "\n";
foreach my $entry (@annos) {
    my @row = ();
    foreach my $t (@terms) {
	push(@row, $entry->{$t});
    }
    my $gref = $entry->{Genes};
    my @genelist = map {$_->{Symbol}} @$gref;
    push(@row, join("|", @genelist));
    print STDERR join("\t", @row), "\n" if $verbose > 1;
    print TSV join("\t", @row), "\n";
}
close(TSV);
print STDERR "Done!\n";
print LOG "Done!\n";

close(LOG);
exit;
__END__
