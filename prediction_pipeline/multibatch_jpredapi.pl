#!/usr/bin/env perl 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                                       #
#   multibatch_jpredapi.pl schedules job submissions to the  REST API of JPred4.        #
#                                                                                       #
#   Copyright (C) 2015                                                                  #
#                                                                                       #
#   ETH Zuerich                                                                         #
#   Institute of Molecular Health Sciences                                              #
#   Fabian EGLI <fabian.egli@biol.ethz.ch)                                              #
#                                                                                       #
#   This program is free software: you can redistribute it and/or modify                #
#   it under the terms of the GNU General Public License as published by                #
#   the Free Software Foundation, either version 3 of the License, or                   #
#   (at your option) any later version.                                                 #
#                                                                                       #
#   This program is distributed in the hope that it will be useful,                     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                       #
#   GNU General Public License for more details.                                        #
#                                                                                       #
#   You should have received a copy of the GNU General Public License                   #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.               #
#                                                                                       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use strict;
use warnings;

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# JPred specific variables. 
my $download_url_base = "http://www.compbio.dundee.ac.uk/jpred4/results/";
my $max_jobs = 1000; # the limit on JPRED REST API per day; 1000 as of March 2015
my $sequlen_max = 800;
my $sequlen_min = 20;
my $checkEvery = 60; # might be anything from 60 to inf
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# set defaults
my $jpredapi = '';
my $verbosity = 0;
my $max_jobs_concurrent = 40; # the limit on JPRED REST API per day; 1000 as of March 2015
my $accnum_list_file = '';
my $email = '';
my $fasta_directory = '';
my $download_directory = 'downloads';
my $download_results = '';
my @download_types;
my $download_file_default = 'tar.gz';
my $logfile = '';
my $logged_job_states = '';
my $print_usage = 0;
my $rename = 0;
my $renamedfolder_name = "renamed_files";
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# parse command line arguments
foreach my $i (0 .. $#ARGV) {
	if ($ARGV[$i] eq '-f') {
		$accnum_list_file = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq '-e') {
		$email = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq '-j') {
		$jpredapi = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq '-i') {
		$fasta_directory = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq '-d') {
		$download_directory = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq '-l') {
		$logfile = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq '-s') {
		$logged_job_states = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq '-c') {
		if ($ARGV[$i+1] > $checkEvery ) {
		    $checkEvery = $ARGV[$i+1];
		}
		else {
		    print "Checking job state should not be done more than once per minute. ceckEvery is set to 60s.\n";
		}
	}
	elsif ($ARGV[$i] eq '-m') {
		$max_jobs_concurrent = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq '-r') {
		$download_results = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq '-h') {
		$print_usage = 1;
	}
	elsif ($ARGV[$i] =~ m/-([v]+)/) {
		$verbosity = length($1);
	}
	elsif ($ARGV[$i] =~ m/-R/) {
		$rename = 1;
		$renamedfolder_name = $ARGV[$i+1] unless ($ARGV[$i+1] =~ m/(^-)|(^$)/);
	}
}

if ($print_usage) {
	print &usage();
	exit();
}

# check command line argument integrity
if ( $jpredapi eq '' ) {
	&print_usage_and_exit("Please provide the jpredapi instance you want to use.\n");
}
unless ( -f $jpredapi ) {
	&print_usage_and_exit("The jpredapi provided does not exist: '$jpredapi'.\n");
}

if ($accnum_list_file eq '') {
	&print_usage_and_exit("Please provide a file containing uniprot accession numbers.\n");
}
unless ( -f $accnum_list_file ) {
	&print_usage_and_exit("The file containing uniprot accession numbers does not exist: '$accnum_list_file'.\n");
}

if ($fasta_directory eq '') {
	&print_usage_and_exit("Please provide a the directory containing the fasta files.\n");
}
unless ( -d $fasta_directory ) {
	&print_usage_and_exit("The fasta directory '$fasta_directory' could not be found.\n$!");
}

if ($email !~ m/^(\S+)@(\S+)\.(\S+)$/) { # A rudimentary email check.
	&print_usage_and_exit("Plese provide a valid email address or none. '$email' is not a valid email adress.");
}

unless ( -d $download_directory ) {
	mkdir($download_directory) or &print_usage_and_exit("Could not create download directory '$download_directory'.\n$!");
}

if ($logfile eq '') {
    $logfile = "$accnum_list_file.log";
}

if ($logged_job_states eq '') {
    $logged_job_states = "$accnum_list_file.job_states.log";
}

if ($download_results eq '') {
    push(@download_types, $download_file_default);
}
else {
    my @download_types_in = split(':', $download_results); # Allowed endings can be found in is_allowed_download_type 
    foreach my $download_type_in (@download_types_in) {
        if ( &is_allowed_download_type($download_type_in) ) {
            push( @download_types, $download_type_in);
        }
        else {
            &print_usage_and_exit("The file suffix '$download_type_in' is not supported.\n");
        }
    }
}

if (! -f $logfile) {
    open my $NL, "> $logfile" or &print_usage_and_exit("Could not create log file $logfile.\n");
    close $NL;
}


my @accnums;
open my $LIST, "< $accnum_list_file" or die "Could not open file contianing accession numbers: '$accnum_list_file'\n$!";
while (my $line = <$LIST>) {
	chomp($line);
	next if ($line eq '');
	push(@accnums, $line);
}
close $LIST;


# defining global variables and setting the defaults for the jobs to be handled by this script
# job status: -1=job failed; 0=job finished; 1=new job; 2=job running
my %job_status = map { $_ => 1} @accnums; # assign all accnums the status new
# all possible jobs get an empty string as an id
my %job_id = map { $_ => ''} @accnums; # assign all accnums the status new
# job_last_checked: 0=never; positive integer=time
my %job_last_checked = map { $_ => 0} @accnums; # assign all accnums the the last checked time of the zero time of Unix time
# job_fail_count: job failed 0 times in the beginning.
my %job_fail_count = map { $_ => 0} @accnums; # assign all jobs a fail count of 0
# job_downloaded: 0=not downloaded; 1=downloaded.
my %job_downloaded = map { $_ => 0} @accnums; # assign all jobs the status not downloaded


###
# inform the user of the script parameters
print "
$0 runs with the following parameters:
    file containing accession numbers: $accnum_list_file
    your email address: $email
    directory containing the files: $fasta_directory
    directory for downloaded predictions: $download_directory
    verbosity $verbosity

";


# parse state and log file
if ( -f $logged_job_states ) {
    &parse_job_states_log($logged_job_states);
}
elsif ( -f $logfile ) {
    &parse_log($logfile);
}
else {
    print "\n -> No log files found, assuming new submission.\n\n";
}


if ($rename) {
	&copy_files;
	&print_statistics if $verbosity;
	exit;
}

# from here on, we want some better handling of programm interruption for consistent log.
$SIG{'INT'} = 'INT_handler';
$SIG{__DIE__}  = 'DIE_handler';

my $remaining_quota_today = &remaining_quota_today();
my $number_of_jobs_to_submit = &number_of_jobs_to_submit();
while (0 < $number_of_jobs_to_submit) {
	print "\nUnsubmitted jobs ...\n" if $verbosity > 2;
	while ( (0 < $remaining_quota_today) and (0 != $number_of_jobs_to_submit) ) {
        print "\nSubmitting jobs ...\n" if $verbosity > 2;
		# Wait until the next submission is allowed
		while ( &number_of_jobs_running() >= $max_jobs_concurrent ) {
			print "Waiting for a job to finish.\n";
			wait_for_one_job_to_finish_and_download();
		}
		
		my $accnum = &next_job_for_submission();
		
		if (defined $accnum) {
            # try submitting the new job
            trysubmit($accnum);
        
            $remaining_quota_today = &remaining_quota_today();
            print &day_time_string_formated_gmt() . ": remaining_quota_today=$remaining_quota_today    accnum=$accnum\n";
		} 
		else {
		    print "All jobs submitted.\n";
		    $number_of_jobs_to_submit = 0;
		}
	}
	
	while ( 0 == $remaining_quota_today and &number_of_jobs_running() > 0) {
		wait_for_one_job_to_finish_and_download();
		$remaining_quota_today = &remaining_quota_today();
	}
	
	if ( $remaining_quota_today == 0 and &number_of_jobs_to_submit() > 0 ) {
		my $sec_to_next_day = &time_left_before_job_quota_re_fill();
		print "Waiting $sec_to_next_day for the new day to start.\n";
		sleep $sec_to_next_day;
		$remaining_quota_today = &remaining_quota_today();
	}
}

wait_for_all_jobs_to_finish_and_download();

print_final_log("EXIT: $!\nEXIT");

exit;




sub trysubmit {
	my ($try_accnum) = @_;
	my $seq_len = &get_seq_length_from_fasta_file($try_accnum);
	if ($job_status{$try_accnum} == 1) {
		if ($seq_len < $sequlen_min) {
			# sequence to long
			$job_status{$try_accnum} = -1;
			print_log("ERROR: $try_accnum: Sequence too short for prediction ($seq_len amino acids).");
		}
		elsif ($seq_len > $sequlen_max) {
			# sequence to long
			$job_status{$try_accnum} = -1;
			print_log("ERROR: $try_accnum: Sequence too long for prediction ($seq_len amino acids).");
		}
		else {
			# submit job
			my ($submit_job_status, $submit_job_job_id, $submit_job_tme, $submit_job_error_message) = &submit_job($try_accnum);
			$job_status{$try_accnum} = $submit_job_status;
			if ($submit_job_status == 2) {
				$job_id{$try_accnum} = $submit_job_job_id;
				$job_last_checked{$try_accnum} = $submit_job_tme;
				print_log("Submitted job for $try_accnum:$submit_job_job_id");
			}
			elsif ($submit_job_status == -1) {
				print_log("ERROR: $try_accnum: job for $try_accnum was rejected: $submit_job_error_message");
				$job_fail_count{$try_accnum}++;
			}
		}
	}
}


sub submit_job {
	my ($accnum) = @_;
	print "submitting job for secondary structure prediction for $accnum\n";
	
	my $fastafile = "$fasta_directory/$accnum.fasta";
	my $logfile = "./log/$accnum.log";
	my $jobname = $accnum;
	
	# Prevent the following error: ERROR: name parameter could only be built from Latin characters, numbers, and '_' symbol. Exiting.
	$jobname =~ s/[^A-Za-z0-9]/\_/g;
	
	my $job_id;
	my $time = 0;
	my $error = undef;
	
	my $status = 1;
	
	# start the job
	print "perl $jpredapi submit file=$fastafile mode=single format=fasta email=$email name=$jobname skipPDB=on", "\n";
	my $cmd_out = `perl $jpredapi submit file=$fastafile mode=single format=fasta email=$email name=$jobname skipPDB=on`;
	# print "cmd_out: $cmd_out\n";
	
	if ($cmd_out =~ m/Created JPred job with jobid: (jp_\S+)/){
		$job_id = $1;
		print "Created job with JobId $job_id.\n";
		$status = 2;
	}
	elsif ($cmd_out =~ m!(ERROR: name parameter could only be built from Latin characters, numbers, and '_' symbol. Exiting.)!){
		$status = -1;
		$error = $1;
	}
	elsif ($cmd_out =~ m!(ERROR: Unrecognised character '(.)' found in sequence\. Only the recognised IUPAC one-letter residue codes are allowed\.This can be caused by selecting the wrong input type for file upload. Please check your setttings and try again.)!){
		$status = -1;
		$error = $1;
		print "Unrecognised character '$2' in input file.\n";
	}
	else {
		print_log("ERROR: unexpected jpredapi output:\n$cmd_out");
	}
	
	$time = time;
	return ($status, $job_id, $time, $error);
}


sub download_results {
	my ($accnum) = @_;
	my $this_job_id = $job_id{$accnum};
	
	my $download_status = '';
	
	foreach my $download_type (@download_types){
        
        my $download_destination = "$download_directory/$this_job_id.$download_type";
		my $download_url = "$download_url_base$this_job_id/$this_job_id.$download_type";
        
        unless ( (-f $download_destination) or ($job_fail_count{$accnum} > 0) ) {
            # do we want the fasta header for the Jnet prediction??? system("head -n 1 $fasta_directory/$accnum.fasta >> $download_file");
            print "curl -s $download_url\n";
            my $jpred_out = `curl -s $download_url`;
        
            if ($jpred_out =~ m/The requested URL is not found\./ ) {
                my $jpred_log = `curl -s $download_url_base$this_job_id/LOG`;
                if ($jpred_log =~ m/--TIMEOUT your job timed out\./ ) {
                    $download_status = 'failed. Job timed out';
                    $job_status{$accnum} = -1;
                    $job_fail_count{$accnum}++;
                }
                elsif ($jpred_log =~ m/(Jpred error:[^\n\r]+)/ ) {
                    $download_status = 'failed.\n$1\n';
                    $job_status{$accnum} = -1;
                    $job_fail_count{$accnum}++;
                }
                else {
                    print "unexpected JPred log:\n$jpred_log\n";
                }
                $download_status = 'fialed';
            }
            else {
                open my $D, "> $download_destination" or die("Could not write the downloaded sequence to the file.");
                print $D $jpred_out;
                close $D;
                $job_downloaded{$accnum} = 1;
                # print "Successfully downloaded $this_job_id.jnet file for $accnum.\n";
                $download_status = 'done';
            }
        }
        else {
            $download_status = 'was already done';
        }
    }
	return $download_status;
}


sub next_job_for_submission {
	# returns the number of jobs waiting for submission
	my @jobs_for_submission = grep {$job_status{$_} == 1} keys %job_status;
	my $jobs_for_submission = shift @jobs_for_submission;
	return $jobs_for_submission;
}


sub wait_for_all_jobs_to_finish_and_download {
	# wait for a job to finish
	while ( (&number_of_jobs_finished() != &number_of_jobs_finished_downloaded()) or (&number_of_jobs_running > 0) ) {
		my $next_check_acc = &get_oldest_running_job_acc();
        
        if (! defined $next_check_acc) {
            $next_check_acc = &get_oldest_finisched_but_not_downloaded_job_acc();
        }
        
		my $time_until_check = $job_last_checked{$next_check_acc} + $checkEvery - time;
		#print "Waiting $time_until_check seconds to check the status of $next_check_acc:$job_id{$next_check_acc}.\n";
		if ($time_until_check > 0) {
			print "Waiting $time_until_check sec before checking if a job finished.\n" if ($verbosity > 0);
			sleep $time_until_check;
		}
		my ($status, $time, $is_success) = &check_job_status($next_check_acc);
		if ($is_success == 0) {
			next;
		}
		# print "status($next_check_acc:$job_id{$next_check_acc}):$status\n";
		if (defined $status) {
			$job_status{$next_check_acc} = $status;
			$job_last_checked{$next_check_acc} = $time;
		
			if ($status == 0) {
				print "downloading $next_check_acc:$job_id{$next_check_acc} ...";
				my $download_status = &download_results($next_check_acc);
				print " $download_status.\n";
			}
		}
	}
	print &number_of_jobs_running() . "\n";
	
}


sub wait_for_one_job_to_finish_and_download {
	# wait for a job to finish
	while ($max_jobs_concurrent <= &number_of_jobs_running()) {
		my $next_check_acc = &get_oldest_running_job_acc();
		my $time_until_check = $job_last_checked{$next_check_acc} + $checkEvery - time;
		#print "Waiting $time_until_check seconds to check the status of $next_check_acc:$job_id{$next_check_acc}.\n";
		if ($time_until_check > 0) {
			print "Waiting $time_until_check sec before checking if a job finished." if ($verbosity > 0);
			sleep $time_until_check;
		}
		my ($status, $time, $is_success) = &check_job_status($next_check_acc);
		if ($is_success == 0) {
			next;
		}
		# print "status($next_check_acc:$job_id{$next_check_acc}):$status\n";
		if (defined $status) {
			$job_status{$next_check_acc} = $status;
			$job_last_checked{$next_check_acc} = $time;
		
			if ($status == 0) {
				print "downloading $next_check_acc:$job_id{$next_check_acc} ...";
				my $download_status = &download_results($next_check_acc);
				print " $download_status.\n";
			}
		}
	}
}


sub check_if_files_downloaded {
	# returns 1 if file exists
	my ($this_acc) = @_;
	my $this_job_id = $job_id{$this_acc};
	# ToDo: check for all files if they are downloaded
	my $download_destination = "$download_directory/$this_job_id.$download_types[0]";
	if (-f $download_destination) {
		#print "Already downloaded: $this_job_id\n";
		$job_status{$this_acc} = 0;
		$job_downloaded{$this_acc} = 1;
		return 1;
	}
	else{
		return 0;
	}
}


sub get_oldest_running_job_acc {
	# returns the job acc that was not checked for the longest time
	my @running_jobs_accs = grep { $job_status{$_} == 2} keys %job_status;
	my @jobs_checked = sort { $job_last_checked{$a} <=> $job_last_checked{$b} } @running_jobs_accs;
	my $acc = shift @jobs_checked;
	return $acc;
}


sub get_oldest_finisched_but_not_downloaded_job_acc {
	# returns the job acc that was not checked for the longest time
	my @running_jobs_accs = grep { $job_status{$_} == 0} keys %job_status;
	my @jobs_checked = sort { $job_last_checked{$a} <=> $job_last_checked{$b} } @running_jobs_accs;
	my $acc;
	while (my $checked = <@jobs_checked>) {
        print "checked: $checked\n";
	    if ($job_downloaded{$checked} == 0) {
	        $acc = $checked;
	        last;
	    }
	}
	return $acc;
}


sub check_job_status {
	my ($job_acc) = @_;
	my $job_id = $job_id{$job_acc};
	my $status;
	my $status_check_is_success;
	
	if ( $job_downloaded{$job_acc} == 1) {
		$status = 0;
	}
	elsif ( &check_if_files_downloaded($job_acc) > 0) {
		$status = 0;
	}
	elsif ( $job_fail_count{$job_acc} > 0) {
		$status = -1;
	}
	else {
		my $cmd_out = `perl $jpredapi status jobid=$job_id checkEvery=once getResults=no`;

		if ($cmd_out =~ m!http://www.compbio.dundee.ac.uk/jpred4/results/$job_id/$job_id.results.html!) {
			$status_check_is_success = 1;
			$status = 0;
		}
		elsif ($cmd_out =~ m/--->	There are currently \d+ jobs due to run before yours\./) {
			$status_check_is_success = 1;
			$status = 2;
		}
		elsif ($cmd_out =~ m/--->\sYour job is next to be submitted\./) {
			$status_check_is_success = 1;
			$status = 2;
		}
		elsif ($cmd_out =~ m/--->\sYour\sjob\s+is\s+\d+%\s+complete\.\.\./) {
			$status_check_is_success = 1;
			$status = 2;
		}
		elsif ($cmd_out =~ m/--->\s+(No job of that ID \($job_id\) was found)\..*/) {
			$status_check_is_success = 1;
			$status = 0;
			print "WARNING: $1:\n$cmd_out\n";
		}
		elsif ($cmd_out =~ m/^\s*\w{3}\s+\w{3}\s+\d{2}\s+\d{2}:\d{2}:\d{2}\s\d{4}\s+--->\s+$/) {
			$status_check_is_success = 0;
			print "WARNING: $1:\n$cmd_out\n";
		}
		else {
			$status_check_is_success = 0;
			print "WARNING: unexpected return value for status check of job '$job_id':\n$cmd_out\n";
		}
	}
	return ($status, time, $status_check_is_success);
}


sub remaining_quota_today {
	my $jpred_quota_query = `perl $jpredapi quota email=$email`;
	my $remaining_quota_today = undef;
	if ( $jpred_quota_query =~ /^You've already submitted (\d+) today \(out of (\d+) jobs per user per day\).$/ ) {
		$remaining_quota_today = $2 - $1;
	}
	elsif ($jpred_quota_query =~ m/You haven't submitted any jobs today yet. If you certain you did submit some jobs - please check spelling of your email./) {
		$remaining_quota_today = $max_jobs;
	}
	else {
		print "ERROR: 'perl $jpredapi quota email=$email' returned unexpected results:\n$jpred_quota_query\n";
	}
	return $remaining_quota_today;
}


sub get_seq_length_from_fasta_file {
	my ($accnum) = @_;
	my $fasta_file = "$fasta_directory/$accnum.fasta";
	
	open my $FILE, "<$fasta_file" or die "Could not open file '$fasta_file'.\n$!";
	
	my $seq_len = 0;
	while (<$FILE>) {
		next if $_ =~ /^>/;
		chomp;
		$seq_len += length;
	}
	return $seq_len;
}


sub INT_handler {
	print "\nCaught a Ctrl-C\n";
    print_final_log("INT: $!\nINT");
    exit();
}


sub DIE_handler {
    my ($string) = @_;
	print "\nCaught a die.\n";
    print "$string\n";
    print_final_log("DIE: $!\nDIE");
    exit();
}


sub print_final_log {
    my($log_string) = @_;
	my @running_jobs_accs = keys %job_status;
	my @running_jobs_acc_id;
	foreach my $running_job (@running_jobs_accs) {
		#print "running_job: $running_job\n";
		my $per_acc_log_str .= '';
		$per_acc_log_str .= "$running_job";
		print "running_job: $running_job\n" unless (defined $job_id{$running_job});
		$per_acc_log_str .= ":$job_id{$running_job}";
		$per_acc_log_str .= ":$job_status{$running_job}";
		$per_acc_log_str .= ":$job_last_checked{$running_job}";
		$per_acc_log_str .= ":$job_downloaded{$running_job}";
		$per_acc_log_str .= ":$job_fail_count{$running_job}";
		push(@running_jobs_acc_id, $per_acc_log_str);
	}
	# print "Running JPred jobs: @running_jobs_acc_id\n"; 
	if (scalar @running_jobs_acc_id > 0) {
		$log_string .= ": job states: " . join(' ', @running_jobs_acc_id);
		print_log($log_string);
		print_job_states_log($log_string);
	}
	else {
		print "No JPred jobs running.\n@running_jobs_accs\n";
		$log_string = "No JPred jobs running.\n";
		print_log($log_string);
		print_job_states_log($log_string);
	}
	
	&print_statistics();

}


sub print_log {
    my($string) = @_;
	open my $LOG, ">> $logfile";
    print $LOG ("$string\n");
    close $LOG;
}


sub print_job_states_log {
    my($string) = @_;
	open my $JS, "> $logged_job_states";
    print $JS ("$string\n");
    close $JS;
}


sub parse_job_states_log {
	my ($logfile) = @_;
	print "Parsing job states log file '$logfile' ...";
	open my $OLD_LOG, "< $logfile" or die "Could not open log file: '$logfile'\n$!";
	while (my $line = <$OLD_LOG>) {
		if ($line =~ m/^(?:INT|DIE|EXIT): job states:\s+(.*)$/) {
			my @logged_running_jobs = split(' ', $1);
			foreach my $logged_running_job (@logged_running_jobs) {
				next if ($logged_running_job =~ /^\s*$/);
				my ($logged_acc, $logged_jobid, $logged_status, $logged_lastchecked, $logged_downloaded, $logged_fail_count) = split(':', $logged_running_job);
				
				print "$logged_acc, $logged_jobid, $logged_status, $logged_lastchecked, $logged_downloaded, $logged_fail_count\n" if $verbosity > 1;
				
				if ( not defined $job_status{$logged_acc} ){
					print "invalid log entry:\n$logged_running_job\n";exit();
				}
				
				if ( $logged_status eq '' ) {
					print "WARNING: Unexpected logged_status in $logged_acc: $logged_status\n";
					my ($chk_status, $chk_time, $chk_is_success) = &check_job_status($logged_acc);
					if ( $chk_is_success ) {
						$job_status{$logged_acc} = $chk_status;
					}
					else {
						# if the job status cn not be found or retrieved from the log just regard it as new and run it again.
						print "WARNING: Unexpected logged_status in $logged_acc: $logged_status\nThis job will be submitted again.\n";
						$job_status{$logged_acc} = 1;
					}
				}
				elsif ($job_status{$logged_acc} == 1) {
					$job_status{$logged_acc} = $logged_status;
					$job_id{$logged_acc} = $logged_jobid;
					$job_last_checked{$logged_acc} = $logged_lastchecked;
					$job_downloaded{$logged_acc} = $logged_downloaded;
					$job_fail_count{$logged_acc} = $logged_fail_count;
				}
				elsif ($job_status{$logged_acc} > $logged_status) {
					$job_status{$logged_acc} = $logged_status;
					$job_id{$logged_acc} = $logged_jobid;
					$job_last_checked{$logged_acc} = $logged_lastchecked;
					$job_downloaded{$logged_acc} = $logged_downloaded;
					$job_fail_count{$logged_acc} = $logged_fail_count;
				}
				else {
					print "invalid log entry:\n$logged_running_job\n" if $verbosity > 0;
				}
			}
		}
		else {
			print "\n The fillowing line is not being parsed:\n$line\n" if $verbosity > 1;
		}
	}
	print " done.\n\n";
	close $OLD_LOG;
}


sub parse_log {
	my ($logfile) = @_;
	print "Parsing old log file '$logfile' ...";
	open my $OLD_LOG, "< $logfile" or die "Could not open log file: '$logfile'\n$!";
	while (my $line = <$OLD_LOG>) {
		if ($line =~ m/^ERROR: ([^\:]+?)\:/){
			if ($line =~ m/^ERROR: unexpected jpredapi output:$/){
				next;
			}
			$job_status{$1} = -1;
		}
		elsif ($line =~ m/^(?:INT|DIE|EXIT): job states:\s+(.*)$/){
			my @logged_running_jobs = split(' ', $1);
			foreach my $logged_running_job (@logged_running_jobs) {
				next if ($logged_running_job =~ /^\s*$/);
				my ($logged_acc, $logged_jobid, $logged_status, $logged_lastchecked, $logged_downloaded, $logged_fail_count) = split(':', $logged_running_job);
				
				if ( not defined $job_status{$logged_acc} ){
					print "invalid log entry:\n$logged_running_job\n";exit();
				}
				
				if ( $logged_status eq '' ){
					print "WARNING: Unexpected logged_status in $logged_acc: $logged_status\n";
					my ($chk_status, $chk_time, $chk_is_success) = &check_job_status($logged_acc);
					if ( $chk_is_success ) {
						$job_status{$logged_acc} = $chk_status;
					}
					else {
						# if the job status cn not be found or retrieved from the log just regard it as new and run it again.
						print "WARNING: Unexpected logged_status in $logged_acc: $logged_status\nThis job will be submitted again.\n";
						$job_status{$logged_acc} = 1;
					}
				}
				elsif ($job_status{$logged_acc} == 1) {
					$job_status{$logged_acc} = $logged_status;
					$job_id{$logged_acc} = $logged_jobid;
					$job_last_checked{$logged_acc} = $logged_lastchecked;
					$job_downloaded{$logged_acc} = $logged_downloaded;
					$job_fail_count{$logged_acc} = $logged_fail_count;
				}
				elsif ($job_status{$logged_acc} > $logged_status) {
					print "T: $logged_status\n";
					$job_status{$logged_acc} = $logged_status;
					$job_id{$logged_acc} = $logged_jobid;
					$job_last_checked{$logged_acc} = $logged_lastchecked;
					$job_downloaded{$logged_acc} = $logged_downloaded;
					$job_fail_count{$logged_acc} = $logged_fail_count;
				}
				else {
					#print "invalid log entry:\n$logged_running_job\n";
					#exit();
				}
			}
		}
		elsif ($line =~ m/^Submitted job for ([^:]+?):(jp_\S+)$/){
			$job_id{$1} = $2;
			if ( &check_if_files_downloaded($1) == 1 ) {
				$job_status{$1} = 0;
				$job_downloaded{$1} = 1;
			}
			else {
				$job_status{$1} = 2;
			}
		}
		else {
			# print $line;
		}
	}
	print " done.\n\n";
	close $OLD_LOG;
}


sub day_time_string_formated_gmt {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time); 
	# print "$sec,$min,$hour,$mday,$mon," . ($year + 1900) . ",$wday,$yday,$isdst\n"
	my @month_name = qw(Jan Feb Mar Apr May Jun Jul Aug Sept Oct Nov Dec);
	my $day_time_string_formated = ($year + 1900) . " $month_name[$mon] $mday, $hour $min' $sec''";
	return $day_time_string_formated;
}


sub time_left_before_job_quota_re_fill {
	my $response = `perl $jpredapi sectonewday`;
	if ( $response =~ /Time left before job quota re-fill: (\d+) \[sec\]/) {
		return $1;
	}
	else {
		print "WARNING: Unexpected return value for refill time check.\n$response\n"; 
	}
}


sub number_of_jobs_to_submit {
	# returns the number of jobs waiting for submission
	my $number_of_jobs_to_submit = grep {$_ == 1} values %job_status;
	return $number_of_jobs_to_submit;
}


sub number_of_jobs_running {
	# returns the number of jobs waiting for submission
	my $number_of_jobs_running = grep {$_ == 2} values %job_status;
	return $number_of_jobs_running;
}


sub number_of_jobs_finished {
	# returns the number of jobs waiting for submission
	my $number_of_jobs_finished = grep {$_ == 0} values %job_status;
	return $number_of_jobs_finished;
}


sub number_of_jobs_finished_downloaded {
	# returns the number of jobs waiting for submission
	my $number_of_jobs_finished_downloaded = grep { $_ == 1 } values %job_downloaded;
	return $number_of_jobs_finished_downloaded;
}


sub number_of_jobs_failed {
	# returns the number of jobs waiting for submission
	my $number_of_jobs_failed = grep { $_ > 0 } values %job_fail_count;
	return $number_of_jobs_failed;
}


sub print_statistics {
	# collect data
	my $job_remaining_count = &number_of_jobs_to_submit();
	my $job_running_count = &number_of_jobs_running();
	my $job_finished_count = &number_of_jobs_finished();
	my $job_downloaded_count = &number_of_jobs_finished_downloaded();
	my $job_fail_count = &number_of_jobs_failed();
	my $jpred_quota = &remaining_quota_today();
	
	my $statistics = "\nStatistics:\n";
	$statistics .= "  jobs:\n";
	$statistics .= "    remaining: $job_remaining_count\n";
	$statistics .= "    running: $job_running_count\n";
	$statistics .= "    finished: $job_finished_count\n";
	$statistics .= "    downloaded: $job_downloaded_count\n";
	$statistics .= "    fialed: $job_fail_count\n";
	$statistics .= "  JPred:\n";
	$statistics .= "    remaining quota: $jpred_quota\n";
	$statistics .= "\n";
	
	print $statistics;
}


sub is_allowed_download_type {
    
    # Readme File as of 13 April 2015
    my $JPred_readme = "
    Filename                                    Description
    --------                                    -----------
    *.align                            PSIBLAST alignment with gaps and redundancy removed in FASTA format
    *.als                              Alscript command file. Used to generate PS/PDF output*.blast.gz                         PSIBLAST output (compressed)
    *.coils.csv                        The output from coils in CSV format
    *.coilseq.lupas_14                 The output from coils using a window length of 14
    *.coilseq.lupas_21                 The output from coils using a window length of 21
    *.coilseq.lupas_28                 The output from coils using a window length of 28
    *.concise                          The prediction in pseudo-CSV format, including the coiled-coil prediction, solvent accessiblity and the sequence alignment
    *.concise.blc                      A BLC file of the prediction and alignment
    *.concise.pdf                      A PDF file of the prediction and alignment
    *.concise.ps                       A PostScript file of the prediction and alignment
    *.fasta                            Input query sequence in FASTA format
    *.full_MSA.fasta                   The full multuple sequence alignment before JPred filters. Gaps/insertions shown. FASTA format.
    *.full_MSA.fasta.noGaps.html       Full multuple sequence alignment before JPred filters. Gaps/insertions not shown. Plain HTML.
    *.full_MSA.fasta.withGaps.html     Full multuple sequence alignment before JPred filters. Gaps/insertions shown. Plain HTML.
    *.hmm                              The HHMer2 profile of the alignment
    *.html                             A HTML file of the prediction and alignment
    *.jalview                          A Jalview annotation file to be read in with the .align file to view the predictions in Jalview
    *.jnet                             The output from Jnet
    *.profile                          PSIBLAST profile
    *.pssm                             PSIBLAST PSSM in a format for Jnet
    *.seq                              Your sequence
    *.simple.html                      The brief HTML output of the query sequence and prediction only
    *.svg.html                         A Jalview generated SVG file with summary of the results.
    ";

    my ($suffix) = @_;
    
    my @job_file_suffixes = qw(align als blast.gz coils.csv coilseq.lupas_14 coilseq.lupas_21 coilseq.lupas_28 concise concise.blc concise.fasta concise.pdf concise.ps e4773326 fasta full_MSA.fasta full_MSA.fasta.noGaps.html full_MSA.fasta.withGaps.html hmm html jalview jnet o4773326 pe4773326 po4773326 profile pssm results.html results_jalview.html seq simple.html svg.html tar.gz);
    
    my %job_file_suffixes = map { $_ => 1} @job_file_suffixes;
    
    my $is_allowed_suffix = 0;
    
    if ( $job_file_suffixes{$suffix} eq 1){
        $is_allowed_suffix = 1;
    }
    
    return $is_allowed_suffix;
}


sub usage {
return "
USAGE:
perl $0 
	
	[description] (Default values; given or explained)
	
    Required command line arguments:
        -j  [path to jpredapi]
        -f  [path to file containing accession numbers]
        -i  [directory containing fasta files]

    Optional command line arguments:
        -e  [your email address]
        -d  [directory to which results should be downloaded] (downloads)
        -l  [log file] (-f input value with suffix '.log')
        -s  [log file for job states] (-f input value with suffix '.job_states.log')
        -c  [equivalent to checkEvery in JPred REST API] (60)
        -m  [maximumal number of concurrent jobs on JPred] (40)
        -r  [':' delimited suffixes for the files that have to be downloaded] (tar.gz)
        -v  [verbosity, you can use -vv for more verbose output]
        -h  print a help message
    
    EXPERIMENTAL: Running the renaming of downloaded files. Currently only works with jnet
    	files.
     	-R	copy and rename the downloaded job files ()

";
}




sub print_usage_and_exit {
    print "\nERROR: $_[0]\n";
    print &usage();
    exit();
}



sub copy_files {
	#"This sub should take all downloaded .jnet files and copy them into a new folder with their names being the fasta input";
	
	print "Copying files ...\n";
	
	if (! -d $renamedfolder_name) {
		mkdir($renamedfolder_name) or die $!;
	}
	
	my @downloaded = grep { $job_downloaded{$_} == 1 } keys %job_downloaded;
	
	foreach my $job (@downloaded) {
		system("cp $download_directory/$job_id{$job}.jnet $renamedfolder_name/$job.jnet");
	}
	
}










