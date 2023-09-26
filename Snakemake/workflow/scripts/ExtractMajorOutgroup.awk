awk '
{
        
    position = 0
    # dna array ACGT
    dna[1] = "A"
    dna[2] = "C"
    dna[3] = "G"
    dna[4] = "T"
    dnaOut = ""
    sumArray[1] = 0
    sumArray[2] = 0
    sumArray[3] = 0
    sumArray[4] = 0
    vsota1 = 0
    vsota2 = 0
    vsota3 = 0
    vsota4 = 0
    for (i=1; i<=NF; i++) {
        # split the line into an array split by the space
        split($i, a, " ")
        # loop through each element splitted by comma and find on which position is number different from 0
        split($i, a, ",")
        for (j=1; j<=4; j++) {
            if (a[j] != 0) {
                position = j
                # increment the sumArray on the position where the number is different from 0
                sumArray[position]++
            }
    }

    # print sum of sumArray[1]
    for (k = 1; k <= sumArray[1]; k++) {
        vsota1 = vsota1 + k
    }
   
    for (k = 1; k <= sumArray[2]; k++) {
        vsota2 = vsota2 + k
    }
    
    for (k = 1; k <= sumArray[3]; k++) {
        vsota3 = vsota3 + k
    }
    
    for (k = 1; k <= sumArray[4]; k++) {
        vsota4 = vsota4 + k
    }
    

    # get highest number from vsota1, vsota2, vsota3, vsota4
    if (vsota1 > vsota2 && vsota1 > vsota3 && vsota1 > vsota4) {
        dnaOut = dna[1]
    }
    if (vsota2 > vsota1 && vsota2 > vsota3 && vsota2 > vsota4) {
        dnaOut = dna[2]
    }
    if (vsota3 > vsota1 && vsota3 > vsota2 && vsota3 > vsota4) {
        dnaOut = dna[3]
    }
    if (vsota4 > vsota1 && vsota4 > vsota2 && vsota4 > vsota3) {
        dnaOut = dna[4]
    }
    # check if there are multiple highest numbers and they are not 0
    if (vsota1 == vsota2 && vsota1 != 0) {
        dnaOut = 0
    }
    if (vsota1 == vsota3 && vsota1 != 0) {
        dnaOut = 0
    }
    if (vsota1 == vsota4 && vsota1 != 0) {
        dnaOut = 0
    }
    if (vsota2 == vsota3 && vsota2 != 0) {
        dnaOut = 0
    }
    if (vsota2 == vsota4 && vsota2 != 0) {
        dnaOut = 0
    }
    if (vsota3 == vsota4 && vsota3 != 0) {
        dnaOut = 0
    }
 
    # reset the sumArray
    sumArray[1] = 0
    sumArray[2] = 0
    sumArray[3] = 0
    sumArray[4] = 0
    }

    print dnaOut

}

' $1
