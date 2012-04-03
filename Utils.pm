# Helper functions used project-wide

package Utils;

use POSIX;

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

1; #Perl requires this for importing a module
