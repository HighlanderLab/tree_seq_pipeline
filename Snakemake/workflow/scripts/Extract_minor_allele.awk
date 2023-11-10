awk '
{
    position = 0
    # dna array ACGT
    dna[1] = "A"
    dna[2] = "C"
    dna[3] = "G"
    dna[4] = "T"
    dnaOut = ""
    min=800
    for (i=1; i<=NF; i++) {
        # split the line into an array split by the comma
        split($i, a, ",")
        # get position of maximum value in array
        for (j=1; j<=length(a); j++) {
            if (a[j] != 0 && a[j] < min) {
                min = a[j]
               	position = j
	        dnaOut = dna[position]  
            }
        }
    }
    print dnaOut
    position = 0
    min = 800
}
' $1
