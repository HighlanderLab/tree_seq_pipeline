awk '
{
    positions[1] = 0
    positions[2] = 0
    positions[3] = 0
    positions[4] = 0
    positions[5] = 0
    positions[6] = 0
    for (i=1; i<=NF; i++) {
        # first postion is position 1, second position is position 2, etc.
        positions[i] = $i
        # check if position 3 contains only one character
        if (i == 3){
            if (length(positions[i]) != 1) {
                # remove the line
                next
            }
        } else if (i == 4){
            # split by comma
            split(positions[i], splitArray, ",")
            # check if every element contains only one character
            for (j=1; j<=length(splitArray); j++) {
                if (length(splitArray[j]) != 1 || splitArray[j] == "*") {
                    # remove the line
                    next
                }
            }
        }
    }

    print positions[1] " " positions[2] " " positions[3] " " positions[4] " " positions[5] " " positions[6]
}

' $1
