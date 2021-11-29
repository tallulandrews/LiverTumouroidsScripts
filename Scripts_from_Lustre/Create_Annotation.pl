use strict;
use warnings;

my $file = "LiverOrganoids-Annotation.csv";

open(my $ifh, $file) or die $!;

my $plate_id = "";
while (<$ifh>) {
	chomp;
	if ($_ =~ /^#(\d+)/) {
		$plate_id = "$1";
		if (length($plate_id) == 3) {
			$plate_id = "0".$plate_id;
		}
	} else {
		my @record = split(/,/);
		my $row = shift(@record);
		for (my $i = 0; $i < scalar(@record); $i++) {
			my $col = $i+1;
			if (length($col) == 1) {
				$col = "0".$i;
			}
			my $well = $row.$col;
			print "$plate_id\_$well\t".$record[$i]."\n";
		}

	}

}
