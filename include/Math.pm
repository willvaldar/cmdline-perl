package Math;
use English;
use StdDefs;
use Nr;
use vars qw( @EXPORT_OK );
use base qw(Exporter);

@EXPORT_OK = qw(
        Apportion
        BinomialCoef
        Factorial
        FactorialLn
        GaussianDeviate
        LinearInterpolate
        LogBase
        Min
        PhiCoef
        PoissonCdf
        PoissonDdf
        PoissonPdf
        PoissonPdfLn
        Round
        Shannon
        SumPairs
        SumPairSeries
        );

use constant BIG    => 1.0e50;
use constant PI     => 3.141592654;
use constant TINY   => 1.0e-20;

sub Apportion($$)
# divides N equitably into n integral parts of size d or d+1.
# It is based on the fact that if
#   n = a + b
# and
#   d = int(N/n)
# then
#   N = ad + b(d+1)
{
    my($total, $nParts) = @ARG;

    if ($total < $nParts)
    {
        Fatal("cannot split N into >N integral parts");
    }

    my $nBigParts = $total % $nParts;
    my $nNormalParts = $nParts - $nBigParts;
    my $normalPartSize = int($total/$nParts);

    return $nNormalParts, $normalPartSize, $nBigParts, $normalPartSize + 1;
}

sub BinomialCoef($$)
{
    return Nr::bico($ARG[0], $ARG[1]);
}

sub ExpCdf($$)
{
    my($lamda, $x) = @ARG;
    my $p = 1 - exp(-$x/$lambda);
    return $p;
}

sub ExpCdfInv($$)
{
    my($lambda, $p) = @ARG;
    my $x = log (1-$p) / -$lambda;
    return $x;
}

sub ExpPdf($$)
{
    my($lamda, $x) = @ARG;
    my $p = $lamda*exp(-$x/$lamda);
    return $p;
}

sub Factorial($)
{
    return Nr::factrl($ARG[0]);
}

sub FactorialLn($)
{
    return Nr::factln($ARG[0]);
}

sub GaussianDeviate($$)
{
    my($mean, $stdDev) = @ARG;

    return Nr::gasdev()*$stdDev + $mean;
}

sub Logit($)
{
	return log($ARG[0]/(1-$ARG[0]));
}

sub LinearInterpolate($$$$)
{
    my($x, $y, $newX, $projectWindow) = @ARG;

    my $maxIndex = scalar(@$x) - 1;

    my $after = $maxIndex;
    for (my $i=0; $i<@$x; $i++)
    {
        if ($newX == $x->[$i])
        {
            return $y->[$i];
        }
        if($newX < $x->[$i])
        {
            $after = $i;
            last;
        }
    }

    if ($after > 0 and $after < $maxIndex)
    {
        my $before = $after - 1;
        my $grad = ($newX - $x->[$before]) / ($x->[$after] - $x->[$before]);
        my $newY = $grad * ($y->[$after] - $y->[$before]) + $y->[$before];
        return $newY;
    }
    elsif (0==$after)
    {
        my $grad = ($y->[$projectWindow-1] - $y->[0]) /
                ($x->[$projectWindow-1] - $x->[0]);
        my $delta = ($x->[0] - $newX)*$grad;
        my $newY = $y->[0] - $delta;
        return $newY;
    }
    elsif ($maxIndex==$after)
    {
        my $grad = ($y->[$maxIndex] - $y->[$maxIndex-$projectWindow+1]) /
                ($x->[$maxIndex] - $x->[$maxIndex-$projectWindow+1]);
        my $delta = ($newX - $x->[$maxIndex])*$grad;
        my $newY = $y->[$maxIndex] + $delta;
        return $newY;
    }
}


sub LogBase ($$)
{
    my($base, $number) = @ARG;

    return log($number) / log($base);
}

sub Min
{
    my $min = shift @ARG;
    for my $n (@ARG)
    {
        $min = $n if $n < $min;
    }
    return $min;
}

sub MultinomDeviate($$)
# WARNING: very computationally inefficient!
# modelled on the R function rmultinom() with n=1
{
    my($size, $aProb) = @ARG;

    my $numBins = scalar @$aProb;

    # initialize cumulative probability data structure
    my $aCumProb = [$aProb->[0]];
    for (my $i=1; $i<$numBins; $i++)
    {
        $aCumProb->[$i] = $aProb->[$i] + $aCumProb->[$i-1];
    }

    my $aCount = Array::New($numBins, 0);
    for (my $i=0; $i<$size; $i++)
    {
        my $r = rand($aCumProb->[(scalar @$aCumProb) -1 ]);
        for (my $p=0; $p<$numBins; $p++)
        {
            if ($r < $aCumProb->[$p])
            {
                $aCount->[$p]++;
                last;
            }
        }
    }
    return $aCount;
}


sub PhiCoef($$$$)
{
    my($p, $n, $o, $u) = @ARG;

    my $numer = ($p*$n - $o*$u);
    my $denom = sqrt( ($p+$o)*($p+$u)*($n+$o)*($n+$u) );

    return (0==$numer) ? 0
                       : $numer / $denom;
}

sub PoissonCdf($$)
{
    my($mean, $value) = @ARG;

    if (0 > $mean or 0 > $value)
    {
        Fatal("cannot have Poisson with mu=$mean, k=$value");
    }
    elsif(0 == $mean)
    {
        return 1;
    }
    elsif (0 == $value)
    {
        return PoissonPdf($mean, 0)
    }
    else
    {
        return Nr::gammq($value, $mean) + PoissonPdf($mean, $value);
    }
}

sub PoissonDdf($$)
{
    my($mean, $value) = @ARG;

    if (0 > $mean or 0 > $value)
    {
        Fatal("cannot have Poisson with mu=$mean, k=$value");
    }
    elsif (0 == $mean)
    {
        return 0==$value ? 1
                         : 0;
    }
    elsif (0 == $value)
    {
        return 1;
    }
    else
    {
        return Nr::gammp($value, $mean);
    }
}

sub PoissonPdf($$)
#_
#_   P(X=k) = e^-m * m^k / k!
#_
{
    my($mean, $k) = @ARG;

    # handle X ~ Poisson(0) seperately to avoid ln(0) error
    if (0 == $mean)
    {
        return (0==$k) ? 1
                       : 0;
    }
    return exp PoissonPdfLn($mean, $k);
}

sub PoissonPdfLn($$)
{
    my($mean, $k) = @ARG;

    # check for invalid parameters
    unless (0 <= $mean and 0 <= $k)
    {
        Fatal("Poisson distribution must have mean and param >= 0."
                ." Mean=$mean and param=$k not allowed");
    }
    # handle X ~ Poisson(0) seperately to avoid ln(0) error
    if (0 == $mean)
    {
        return 0 if 0==$k;
        Fatal("cannot take ln(P) when P=0");
    }

    my $lnP = -$mean + $k*log($mean) - FactorialLn($k);

    return $lnP;
}

sub Round
{
    my($x, $numPlaces) = @ARG;
    $numPlaces = 0 unless defined $numPlaces;

    return sprintf "%.${numPlaces}f ",$x;
}

sub Shannon($)
{
    my $aFractionalFreqs = shift;

    my $sum = 0;
    for my $pi (@$aFractionalFreqs)
    {
        if ($pi > 0)
        {
            $sum += $pi * log $pi;
        }
    }
    return (-$sum/log(2));
}

sub SumPairSeries($)
{
    my $n = shift;

    return $n*($n-1)/2;
}

sub SumPairs($)
{
    my $n = shift;

    warn "deprecated method SumPairs\n";

    return SumPairsSeries($n);
}

true;
