#!/usr/bin/perl 

use strict;
use warnings;
use diagnostics;
use Getopt::Long;

my $N = 0;
my $np = 0;
my $ref_levels = 0;
my $dry = 0;
my $flamegraph = 0;
my $degree = 1;
my $statsfile = 'data_mms.csv';
my $csolver = 'asm';

GetOptions("N=i" => \$N,
"np=i" => \$np,
"ref=i" => \$ref_levels,
"csolver=s" => \$csolver,
"p=i" => \$degree,
"statsfile=s" => \$statsfile,
"flamegraph" => \$flamegraph,
"dry" => \$dry) or die("Error in command line arguments\n");

my $jobid = sprintf("mms_%s_%s_%s_%s_%s", $N, $np, $ref_levels, $csolver, $degree);
my $filename = sprintf("%s.job", $jobid);

my $log_base = sprintf("%s_", $jobid);
my $log_number = 0;
my $logname = $log_base . sprintf("%#.2o", $log_number);
my $logfilename = $logname . '.log';

# while (-e $logfilename) {
#     $log_number = $log_number + 1;
#     $logname = $log_base . sprintf("%#.2o", $log_number);
# }

my $str = <<END;
#!/bin/bash
#SBATCH -N $N
#SBATCH -A sosu
#SBATCH -t 60
#SBATCH --output=$logfilename

/usr/bin/time -v srun -N $N -n $np python mms_scaling.py \\
-ref_levels=$ref_levels \\
-stats_file "$statsfile" \\
-csolver "$csolver" \\
-degree $degree \\
END

if ($flamegraph) {
    $str = $str . "-log_view :flamegraph_data_$jobid:ascii_flamegraph"
} else {
    $str = $str . "-log_view"
}

$str = $str . "\n";

if ($dry == 1) {
    print($str)
}
else {
    open(FH, '>', $filename) or die $!;
    print FH $str;
    close(FH);

    system("rm $logfilename");
    system(sprintf("sbatch %s.job", $jobid));
}
