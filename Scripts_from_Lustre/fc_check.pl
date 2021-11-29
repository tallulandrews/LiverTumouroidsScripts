my @files1 = glob("BAMS/*.bam");
my @files2 = glob("Counts/*counts");
my %ids = ();
foreach my $file (@files2) {
	if ($file =~/(_\d_\d+)\./) {
		$ids{$1} = 1;
	} else {
		die "$file does not match\n";
	}
}

my $count=0;
foreach my $file (@files1) {
	$count++;
	if ($file =~/(_\d_\d+)\./) {
		if (!exists($ids{$1})) {
			print $1."\t".$count."\n";
		}
	} else {
		die "$file does not match\n";
	}
}
