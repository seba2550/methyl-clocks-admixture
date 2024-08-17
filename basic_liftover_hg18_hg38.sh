# Use this to liftover coordinates from hg18 to hg38. Because I have extra columns in the horvath clock file (and it is "technically" no longer in bed format), I have to add the bedPlus parameter and specify that only the first three columns
should be used for liftover.
./liftover -bedPlus=3 Horvath_clock_coefficients.bed hg18ToHg38.over.chain.gz Horvath_clock_coefficients_hg38.bed horvath_clock_hg18_unlifted.bed
