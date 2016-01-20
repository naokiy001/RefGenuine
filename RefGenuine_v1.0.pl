# **** Software tool name: RefGenuine ****
# Platform: Perl
# This tool was developed for finding the best reference genes in gene expression analyses.
# Users have to prepare Ct values of potential reference genes across multiple target samples,
# and the tool incorporates the Ct values to calculate an index of gene expression stability
# by the method of Silver et al. (2006).

# The index will be calculated according to the following initial setting parameter:
#  How many gene do make gene combination? Up to how many genes?
# In current, this tool accepts upto 10 gene combinations at max. Please edit the number of 
# genes in the equation below ($combi). The default value $combi is set as "4". The value has
# to be set no more than your candidate reference genes tested.

$combi = 4 ; # Number of maximum reference genes

# **** File format definition ****
# Input file has to be given in a tab-delimitated text file (please see the example file "singlect.txt")
#   file name: singlect.txt
#   text encoding: SJIS code text
#   The 1st row: header (no processed)
#   The 2nd row later : data content
#   The 1st colurm : Gene identifier (user-defined)
#   The 2nd colourm later : Ct value in each sample
# The input file must be put at the same directory with this tool
#
# Two output files will be created after finishing run.
# 1) AveMeanSDdeltaCt.txt: This represents the gene expression stability index.
# 2) Calculation: This represents meanSD of deltaCt and standard deviation of meanSD of deltaCt (intermediate data).
#
# These files are both tab-deliminated.
# In the file AveMeanSDdeltaCt, the best reference gene shows the minimum value of the index.
###################################################################################################
open (IN, "singlect.txt") or die;
while ($D20 = <IN>) {
	chomp $D20 ;
	$skip++;
	if ($skip == 1) {
		@D20count = split (/\t/, $D20) ;
		$sn = @D20count - 1 ;
		next ;
	}
	$a++ ;
	@D20count = split (/\t/, $D20) ;
	for ($x = 0; $x <= $sn ; $x++) {
		if ($x == 0) {
			$GeneID{$a} = $D20count[$x] ;
		} else {
			$Ct{$a}{$x} = $D20count[$x] ;
		}
	}
}
@D20count = () ;
close (IN) ;
$b = $a ;

if ($combi < 2) {
	goto calc ;
}
for ($x = 1; $x <= $a; $x++) {
	for ($y = ($x + 1) ; $y <= $a; $y++) {
		$b++ ;
		$GeneID{$b} = "$GeneID{$x}_$GeneID{$y}" ;
		for ($z = 1; $z <= $sn ; $z++) {
			$cCt = sqrt($Ct{$y}{$z} * $Ct{$x}{$z}) ;
			$Ct{$b}{$z} = $cCt ;
		}
	}
}

if ($combi < 3) {
	goto calc ;
}
for ($x = 1; $x <= $a; $x++) {
	for ($y = ($x + 1) ; $y <= $a; $y++) {
		for ($z = ($y + 1) ; $z <= $a; $z++) {
			$b++ ;
			$GeneID{$b} = "$GeneID{$x}_$GeneID{$y}_$GeneID{$z}" ;
			for ($i = 1; $i <= $sn; $i++) {
				$cCt = ($Ct{$y}{$i} * $Ct{$x}{$i} * $Ct{$z}{$i}) ** (1/3) ;
				$Ct{$b}{$i} = $cCt ;
			}
		}
	}
}

if ($combi < 4) {
	goto calc ;
}
for ($x = 1; $x <= $a; $x++) {
	for ($y = ($x + 1) ; $y <= $a; $y++) {
		for ($z = ($y + 1) ; $z <= $a; $z++) {
			for ($zz = ($z + 1) ; $zz <= $a; $zz++) {
				$b++ ;
				$GeneID{$b} =  "$GeneID{$x}_$GeneID{$y}_$GeneID{$z}_$GeneID{$zz}" ;
				for ($i = 1; $i <= $sn; $i++) {
					$cCt = ($Ct{$y}{$i} * $Ct{$x}{$i} * $Ct{$z}{$i} * $Ct{$zz}{$i}) ** (1/4) ;
				$Ct{$b}{$i} = $cCt ;
				}
			}
		}
	}
}

if ($combi < 5) {
	goto calc ;
}
for ($x = 1; $x <= $a; $x++) {
	for ($y = ($x + 1) ; $y <= $a; $y++) {
		for ($z = ($y + 1) ; $z <= $a; $z++) {
			for ($xz = ($z + 1) ; $xz <= $a; $xz++) {
				for ($zz = ($xz + 1) ; $zz <= $a; $zz++) {
					$b++ ;
					$GeneID{$b} =  "$GeneID{$x}_$GeneID{$y}_$GeneID{$z}_$GeneID{$xz}_$GeneID{$zz}" ;
					for ($i = 1; $i <= $sn; $i++) {
						$cCt = ($Ct{$y}{$i} * $Ct{$x}{$i} * $Ct{$z}{$i} * $Ct{$xz}{$i} * $Ct{$zz}{$i}) ** (1/5) ;
						$Ct{$b}{$i} = $cCt ;
					}
				}
			}
		}
	}
}

if ($combi < 6) {
	goto calc ;
}
for ($x = 1; $x <= $a; $x++) {
	for ($y = ($x + 1) ; $y <= $a; $y++) {
		for ($z = ($y + 1) ; $z <= $a; $z++) {
			for ($xz = ($z + 1) ; $xz <= $a; $xz++) {
				for ($yz = ($xz + 1) ; $yz <= $a; $yz++) {
					for ($zz = ($yz + 1) ; $zz <= $a; $zz++) {
						$b++ ;
						$GeneID{$b} =  "$GeneID{$x}_$GeneID{$y}_$GeneID{$z}_$GeneID{$xz}_$GeneID{$yz}_$GeneID{$zz}" ;
						for ($i = 1; $i <= $sn; $i++) {
							$cCt = ($Ct{$y}{$i} * $Ct{$x}{$i} * $Ct{$z}{$i} * $Ct{$xz}{$i} * $Ct{$yz}{$i} * $Ct{$zz}{$i}) ** (1/6) ;
							$Ct{$b}{$i} = $cCt ;
						}
					}
				}
			}
		}
	}
}

if ($combi < 7) {
	goto calc ;
}
for ($x = 1; $x <= $a; $x++) {
	for ($y = ($x + 1) ; $y <= $a; $y++) {
		for ($z = ($y + 1) ; $z <= $a; $z++) {
			for ($xz = ($z + 1) ; $xz <= $a; $xz++) {
				for ($yz = ($xz + 1) ; $yz <= $a; $yz++) {
					for ($vz = ($yz + 1) ; $vz <= $a; $vz++) {
						for ($zz = ($vz + 1) ; $zz <= $a; $zz++) {
							$b++ ;
							$GeneID{$b} =  "$GeneID{$x}_$GeneID{$y}_$GeneID{$z}_$GeneID{$xz}_$GeneID{$yz}_$GeneID{$vz}_$GeneID{$zz}" ;
							for ($i = 1; $i <= $sn; $i++) {
								$cCt = ($Ct{$y}{$i} * $Ct{$x}{$i} * $Ct{$z}{$i} * $Ct{$xz}{$i} * $Ct{$yz}{$i} * $Ct{$vz}{$i} * $Ct{$zz}{$i}) ** (1/7) ;
								$Ct{$b}{$i} = $cCt ;
							}
						}
					}
				}
			}
		}
	}
}

if ($combi < 8) {
	goto calc ;
}
for ($x = 1; $x <= $a; $x++) {
	for ($y = ($x + 1) ; $y <= $a; $y++) {
		for ($z = ($y + 1) ; $z <= $a; $z++) {
			for ($xz = ($z + 1) ; $xz <= $a; $xz++) {
				for ($yz = ($xz + 1) ; $yz <= $a; $yz++) {
					for ($vz = ($yz + 1) ; $vz <= $a; $vz++) {
						for ($uz = ($vz + 1) ; $uz <= $a; $uz++) {
							for ($zz = ($uz + 1) ; $zz <= $a; $zz++) {
								$b++ ;
								$GeneID{$b} =  "$GeneID{$x}_$GeneID{$y}_$GeneID{$z}_$GeneID{$xz}_$GeneID{$yz}_$GeneID{$vz}_$GeneID{$uz}_$GeneID{$zz}" ;
								for ($i = 1; $i <= $sn; $i++) {
									$cCt = ($Ct{$y}{$i} * $Ct{$x}{$i} * $Ct{$z}{$i} * $Ct{$xz}{$i} * $Ct{$yz}{$i} * $Ct{$vz}{$i} * $Ct{$uz}{$i} * $Ct{$zz}{$i}) ** (1/8) ;
									$Ct{$b}{$i} = $cCt ;
								}
							}
						}
					}
				}
			}
		}
	}
}

if ($combi < 9) {
	goto calc ;
}
for ($x = 1; $x <= $a; $x++) {
	for ($y = ($x + 1) ; $y <= $a; $y++) {
		for ($z = ($y + 1) ; $z <= $a; $z++) {
			for ($xz = ($z + 1) ; $xz <= $a; $xz++) {
				for ($yz = ($xz + 1) ; $yz <= $a; $yz++) {
					for ($vz = ($yz + 1) ; $vz <= $a; $vz++) {
						for ($uz = ($vz + 1) ; $uz <= $a; $uz++) {
							for ($tz = ($uz + 1) ; $tz <= $a; $tz++) {
								for ($zz = ($tz + 1) ; $zz <= $a; $zz++) {
									$b++ ;
									$GeneID{$b} =  "$GeneID{$x}_$GeneID{$y}_$GeneID{$z}_$GeneID{$xz}_$GeneID{$yz}_$GeneID{$vz}_$GeneID{$uz}_$GeneID{$tz}_$GeneID{$zz}" ;
									for ($i = 1; $i <= $sn; $i++) {
										$cCt = ($Ct{$y}{$i} * $Ct{$x}{$i} * $Ct{$z}{$i} * $Ct{$xz}{$i} * $Ct{$yz}{$i} * $Ct{$vz}{$i} * $Ct{$uz}{$i} * $Ct{$tz}{$i} * $Ct{$zz}{$i}) ** (1/9) ;
										$Ct{$b}{$i} = $cCt ;
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

if ($combi < 10) {
	goto calc ;
}
for ($x = 1; $x <= $a; $x++) {
	for ($y = ($x + 1) ; $y <= $a; $y++) {
		for ($z = ($y + 1) ; $z <= $a; $z++) {
			for ($xz = ($z + 1) ; $xz <= $a; $xz++) {
				for ($yz = ($xz + 1) ; $yz <= $a; $yz++) {
					for ($vz = ($yz + 1) ; $vz <= $a; $vz++) {
						for ($uz = ($vz + 1) ; $uz <= $a; $uz++) {
							for ($tz = ($uz + 1) ; $tz <= $a; $tz++) {
								for ($sz = ($tz + 1) ; $sz <= $a; $sz++) {
									for ($zz = ($sz + 1) ; $zz <= $a; $zz++) {
										$b++ ;
										$GeneID{$b} =  "$GeneID{$x}_$GeneID{$y}_$GeneID{$z}_$GeneID{$xz}_$GeneID{$yz}_$GeneID{$vz}_$GeneID{$uz}_$GeneID{$tz}_$GeneID{$sz}_$GeneID{$zz}" ;
										for ($i = 1; $i <= $sn; $i++) {
											$cCt = ($Ct{$y}{$i} * $Ct{$x}{$i} * $Ct{$z}{$i} * $Ct{$xz}{$i} * $Ct{$yz}{$i} * $Ct{$vz}{$i} * $Ct{$uz}{$i} * $Ct{$tz}{$i} * $Ct{$sz}{$i} * $Ct{$zz}{$i}) ** (1/10) ;
											$Ct{$b}{$i} = $cCt ;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

calc : {
open (OUT1, ">Calculation.txt") ; # Gene1,Gene2, meanSD of deltaCt, Stdev of meanSD deltaCt
open (OUT2, ">AveMeanSDdeltaCt.txt") ;
for ($x = 1; $x <= $b; $x++) {
	$sumsd = $cc = 0 ;
	for ($y = 1; $y <= $b; $y++) {
		$sumdCt = $temp = $sd = 0;
		if ($x == $y) {
			next ;
		} else {
			for ($z = 1; $z <= $sn; $z++) {
				$dCt{$z} = $Ct{$y}{$z} - $Ct{$x}{$z} ;
				$sumdCt = $sumdCt + $dCt{$z} ;
			}
			$avsumdCt = $sumdCt / $sn ;
			for ($z = 1; $z <= $sn; $z++) {
				$temp = $temp + ($dCt{$z} - $avsumdCt) * ($dCt{$z} - $avsumdCt) ;
			}
			$sd = sqrt($temp/$sn) ;
		}
		print OUT1 "$GeneID{$x}\t$GeneID{$y}\t$avsumdCt\t$sd\n" ;
		$sumsd = $sumsd + $sd ;
		$cc++ ;
	}
	$avsumsd = ($sumsd / $cc)  ; 
	print OUT2 "$GeneID{$x}\t$avsumsd\n" ;
}
close (OUT1) ;
close (OUT2) ;
}
