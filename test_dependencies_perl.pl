eval {
    require IO::Uncompress::Gunzip;
    IO::Uncompress::Gunzip->import( qw/gunzip $GunzipError/ ) ;
};
if ($@) {
    die "Error: Perl Module IO::Uncompress::Gunzip not installed";
}