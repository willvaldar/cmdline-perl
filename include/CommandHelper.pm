package VUtils::CommandHelper;

use Exporter;
use FileHandle;
use English;
use Carp;

#_
#_ Pragmata
#_

use strict;
use vars qw( @EXPORT_OK );

#_
#_ Inheritance
#_

use base qw( Exporter );

@EXPORT_OK = qw(
        GetRegExpFromFile
        UserConfirm
    );


sub GetRegExpFromFile ($)
{
    my($regExpFileName) = @ARG;
    
    my $aRegExp = [];
    my $line;
    my $fhre = new FileHandle($regExpFileName, "r")
    or confess "Exception: cannot open file $regExpFileName for reading: $ERRNO\n";
    while($line = $fhre->getline())
    {
        if($line =~ /^\s*(.+)\s*$/)
        {
            push @$aRegExp, $1;
        }
    }
    $fhre->close();
    return $aRegExp;
}

sub UserConfirm ($)
{
    my ($message) = @ARG;
    print $message;
    return (<STDIN> =~ m/^[yY]/);
}
