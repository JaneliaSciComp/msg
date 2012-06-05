# Helper functions used project-wide

package Utils;

use POSIX;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;

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
