#!/usr/bin/perl


open(OUT,">do_all_new.sh");
open(FH,"<get_proc_elines_MUSE.csv");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	my @data=split(",",$line);
	$name=$data[0];
	$red=$data[14];
	$check_file="/disk-d/manga/carlos/ccf_products/MUSE_ccf_200620/ccf_dataproducts/".$name.".Ha_NII.ccf.vel.fits";
	if (!(-e $check_file)) {
	    $call="tsp ./ccf_IFS.py  /nfs/ofelia/disk-b/sanchez/MUSE/AMUSING/analysis/all_files_SSP/".$name."/GAS.".$name.".cube.fits.gz ".$red." 2.6 ".$name.".Ha_NII 0 /disk-d/manga/carlos/ccf_products/MUSE_ccf_200620/ccf_dataproducts/";
	    #	    $call="tsp ./ccf_IFS.py /nfs/ofelia/disk-b/sanchez/MUSE/AMUSING/analysis/all_files_SSP/".$name."/GAS.".$name.".cube.fits.gz ".$red." /disk-d/manga/carlos/ccf_products/MUSE/ccf_dataproducts/".$name.".Ha_NII. 0";
	    print OUT "$call\n";
	} else {
	    print "$name already ccf analyzed\n";
	}
    }
}
close(FH);
close(OUT);

exit;
