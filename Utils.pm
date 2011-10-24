# Helper functions used project-wide

package Utils;

sub strip {
    # Remove leading and trailing whitespace
    my ($val) = @_;
    $val =~ s/^\s+//;
    $val =~ s/\s+$//;
    return $val;
}


1; #Perl requires this for importing a module
