package StdDefs;
use strict;
use vars qw( @ISA @EXPORT );
use Exporter;

@ISA = qw( Exporter );
@EXPORT = qw( true false null );

use constant true => 1;
use constant false => 0;
use constant null => '';

