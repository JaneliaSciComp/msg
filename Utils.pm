# Helper functions used project-wide

package Utils;

use POSIX;
eval {
    require IO::Uncompress::Gunzip;
    IO::Uncompress::Gunzip->import( qw/gunzip $GunzipError/ ) ;
};
if ($@) {
    die "Error: Perl Module IO::Uncompress::Gunzip not installed";
}

sub strip {
    # Remove leading and trailing whitespace
    my ($val) = @_;
    $val =~ s/^\s+//;
    $val =~ s/\s+$//;
    return $val;
}

sub system_call {
    print "\nstarted ".POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime);
    print "  @_\n" ;
    system("@_") == 0 or die "Error in @_: $?" ;
    print 'ended '.POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime);
}

sub test_dependencies {
    # Make sure all required dependencies are installed
    my $last_path = getcwd();
    chdir('msg') or die "$!";
    system_call("chmod 755 test_dependencies.sh");
    system_call("test_dependencies.sh");
    chdir($last_path) or die "$!";
}

sub parse_config {
    #Read in msg.cfg or other specified file and update where needed
    
    my ($cfg_path, $default_params) = @_;
    my %default_params = %$default_params;

    die "ERROR: Can't locate $cfg_path.\n" unless (-e $cfg_path);
    open (IN, $cfg_path) || die "ERROR: Can't open $cfg_path: $!\n";
    while (<IN>) { chomp $_;
        next if ($_ =~ /^\#/);
        next unless ($_);
        my ($key,$val) = split(/=/,$_,2);
        $default_params{strip($key)} = strip($val);
    } close IN;
    
    ### Configure some parameters ###
    if (defined $default_params{'chroms'}) {
        $default_params{'chroms2plot'} = $default_params{'chroms'} unless (defined $default_params{'chroms2plot'});
    }
    my $update_nthreads = $default_params{'threads'} if (defined $default_params{'threads'}); ## Number of qsub slots when running pe option
    #add space after qsub options so we can insert into commands, add thread/slot count to -pe option
    if (defined $default_params{'addl_qsub_option_for_exclusive_node'} && $default_params{'addl_qsub_option_for_exclusive_node'}) {
        #example: go from user msg.cfg entered "-l excl=true" to "-l excl=true "
        $default_params{'addl_qsub_option_for_exclusive_node'} = $default_params{'addl_qsub_option_for_exclusive_node'}.' ';
    }
    else {
        $default_params{'addl_qsub_option_for_exclusive_node'} = '';
    }
    if (defined $default_params{'addl_qsub_option_for_pe'} && $default_params{'addl_qsub_option_for_pe'}) {
        #example: go from user msg.cfg entered "-pe batch" to "-pe batch 8 "
        $default_params{'addl_qsub_option_for_pe'} = $default_params{'addl_qsub_option_for_pe'}." $update_nthreads ";
    }
    else {
        $default_params{'addl_qsub_option_for_pe'} = '';
    }
    return \%default_params;
}

sub validate_config {
    #Check that params were entered correctly and nothing essential is missing
    my ($params, @required_file_paths) = @_;
    ### check if all files exist
    foreach my $param (@required_file_paths) {
        if (exists $params->{$param}) {
            die "Exiting from msgCluster: Missing file $params->{$param}.\n" unless (-e $params->{$param});
        }
    }
    print "\nParameters:\n\n";
    ### double check if the minimum exist
    foreach my $key (sort keys %$params) {
        die "ERROR (msgCluster): undefined parameter ($key) in config file.\n" unless ($params->{$key} ne 'NULL');
        print "$key:\t$params->{$key}\n" ;
    }
    print "\n" ;
}

sub readFasta {
    #Return a dictionary: length or sequence of each read by reference id.  Supports gzipped 
    #fasta files or regular.
    my ($file, $store_count) = @_;
    
    my $file_handle;
    my %reads;
    my ($read,$seq);
    
    if ($file =~ /\.gz$/ || $file =~ /\.gzip$/) {        
        $file_handle = new IO::Uncompress::Gunzip $file
            or die "gunzip failed: $GunzipError\n";    
    } else {
        open $file_handle, "<", $file || die "ERROR : Can't open $file: $!\n";
    }

    while (<$file_handle>) { 
        chomp $_;
        if ($_ =~ /^>(\S+)/) {
            #If we don't have a previous read to save then something went wrong, check first
            if ($seq && !$read) {die "Invalid FASTA file: @_";}
            if ($seq) {
                if ($store_count) {
                    $reads{$read} = length($seq)
                } else {
                    $reads{$read} = $seq
                }            
            }
            $read = $1;
            $seq = '';
        } else { $seq .= $_; }
    } 
    close $file_handle;
    
    if ($read) {
        if ($store_count) {
            $reads{$read} = length($seq)
        } else {
            $reads{$read} = $seq
        }
    }    
    return %reads;  
}

1; #Perl requires this for importing a module
