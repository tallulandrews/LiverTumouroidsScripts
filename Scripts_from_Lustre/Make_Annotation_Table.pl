use strict;
use warnings;


if (@ARGV < 1) {die "requires at least one headerfile\n";}
my %WellID2CellID = ();

foreach my $file (@ARGV) {
	#Extract cell ID
	my $cellid = "";
	if ($file =~ /_(\d_\d+)\./) {
		$cellid = $1;
	} else {
		die "$file does not match\n";
	}
	open (my $ifh, $file) or die $!;
	while (<$ifh>) {
		if ($_ =~ /^\@RG/) {
			# Match the plate-well ID
			my $wellid = "";
			if ($_ =~ /SM:SCGC--(\w+)/) {
				$wellid = $1;
			} else {
				die "$_ does not match";
			}
			$WellID2CellID{$wellid} = $cellid;
			last;
		} else {
			next;
		}
	} close($ifh);
}

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
				$col = "0".($i+1);
			}
			my $well = $row.$col;
			my $well_id = "$plate_id\_$well";
			print $WellID2CellID{$well_id}."\t$well_id\t".$record[$i]."\n";
		}

	}

}
