use strict;
use warnings;

my @files1 = glob("/warehouse/team218_wh01/MH/LBHsLiverData/21698/Lane8/*");
my @files2 = glob("/warehouse/team218_wh01/MH/LBHsLiverData/21698/Lane7/*");
my @files3 = glob("/warehouse/team218_wh01/MH/LBHsLiverData/21698/Lane6/*");

my %cell2files = ();
for my $file (@files1, @files2, @files3) {
	if ($file =~ /21698_\d_(\d+).cram/) {
		$cell2files{$1}->{$file} = 1;
	} else {
		die "$file does not match!";
	}
}

foreach my $cell (sort(keys(%cell2files))) {
	print $cell."\t".scalar(keys(%{$cell2files{$cell}}))."\n";
}
