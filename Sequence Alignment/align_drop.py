import math


def main():
    score_matrix, letters, thresholdX, bandwidth, h, s, seq1, seq2, seq_names, outputfile_name = readfile()

    # print inputs
    print("Sequence 1: " + seq1)
    print("Sequence 2: " + seq2)
    print("threshold: % d, bandwidth: % d, starting gap penalty: % d, gap extension penalty: % d" %(thresholdX, bandwidth, h, s))
    print()

    # initialization
    col = len(seq1)
    row = len(seq2)
    entries = row*col
    h = -h
    s = -s
    V = [[0 for i in range(col)] for j in range(row)]
    consistency = [[0 for i in range(col)] for j in range(row)]
    E = [[0 for i in range(col)] for j in range(row)]
    F = [[0 for i in range(col)] for j in range(row)]

    maximum = 0
    m_row = 0
    m_col = 0

    V[0][0] = 0
    consistency[0][0] = 1
    entries = 1

    d = row + col - 1
    # for each diagonal
    for diagonal in range(1, d):
        has_consistent = False
        counter = 0
        if diagonal < row:
            # start diagonal from V[diagonal][0]
            i = diagonal
            j = 0
            # to process diagonal from V[i][j]
            for x in range(col):
                # ---------process V[i][j]----------------
                a_con, b_con, c_con = 0, 0, 0
                if i - 1 >= 0:
                    a_con = consistency[i - 1][j]
                if j - 1 >= 0:
                    b_con = consistency[i][j - 1]
                if i - 1 >= 0 and j - 1 >= 0:
                    c_con = consistency[i - 1][j - 1]
                if max(a_con, b_con, c_con) == 0:
                    # inconsistent entry
                    i -= 1
                    j += 1
                    if i < 0 or j > col - 1:
                        break
                    continue

                E[i][j], F[i][j], aligned = -math.inf, -math.inf, -math.inf
                if i - 1 >= 0:
                    E[i][j] = max(E[i][j-1]-s, V[i][j-1]-h-s)
                    entries += 1
                if j - 1 >= 0:
                    F[i][j] = max(F[i-1][j]-s, V[i-1][j]-h-s)
                    entries += 1
                if i - 1 >= 0 and j - 1 >= 0:
                    aligned = V[i - 1][j - 1] + align(seq1[j], seq2[i], letters, score_matrix)
                V[i][j] = max(aligned, F[i][j], E[i][j])
                entries += 1

                if V[i][j] >= maximum - thresholdX:
                    consistency[i][j] = 1
                    has_consistent = True
                if V[i][j] > maximum:
                    maximum = V[i][j]
                    m_row = i
                    m_col = j
                # -------------------------------------
                i -= 1
                j += 1
                if i < 0 or j > col - 1:
                    break
        else:
            # start diagonal from V[row-1][diagonal-row+1]
            i = row-1
            j = diagonal-row+1
            # to process diagonal from V[i][j]
            for x in range(col):
                # ---------process V[i][j]----------------
                a_con, b_con, c_con = 0, 0, 0
                if i - 1 >= 0:
                    a_con = consistency[i - 1][j]
                if j - 1 >= 0:
                    b_con = consistency[i][j - 1]
                if i - 1 >= 0 and j - 1 >= 0:
                    c_con = consistency[i - 1][j - 1]
                if max(a_con, b_con, c_con) == 0:
                    # inconsistent entry
                    i -= 1
                    j += 1
                    if i < 0 or j > col - 1:
                        break
                    continue

                E[i][j], F[i][j], aligned = -math.inf, -math.inf, -math.inf
                if i - 1 >= 0:
                    E[i][j] = max(E[i][j - 1] - s, V[i][j - 1] - h - s)
                    entries += 1
                if j - 1 >= 0:
                    F[i][j] = max(F[i - 1][j] - s, V[i - 1][j] - h - s)
                    entries += 1
                if i - 1 >= 0 and j - 1 >= 0:
                    aligned = V[i - 1][j - 1] + align(seq1[j], seq2[i], letters, score_matrix)
                V[i][j] = max(aligned, F[i][j], E[i][j])
                entries += 1

                if V[i][j] >= maximum - thresholdX:
                    consistency[i][j] = 1
                    has_consistent = True
                if V[i][j] > maximum:
                    maximum = V[i][j]
                    m_row = i
                    m_col = j
                # -------------------------------------
                i -= 1
                j += 1
                if i < 0 or j > col-1:
                    break

        # if cannot generate anymore entries in  next consecutive diagonals
        if has_consistent:
            counter = 0
        else:
            counter += 1
            if counter > 1:
                break

    # print matrix V
    print("V: ", end=" ")
    for i in range(len(seq1)):
        print(seq1[i], end="  ")
    print()
    for i in range(len(seq2)):
        print(seq2[i], V[i])
    print()

    # print matrix E
    print("E: ", end=" ")
    for i in range(len(seq1)):
        print(seq1[i], end="  ")
    print()
    for i in range(len(seq2)):
        print(seq2[i], E[i])
    print()

    # print matrix F
    print("F: ", end=" ")
    for i in range(len(seq1)):
        print(seq1[i], end="  ")
    print()
    for i in range(len(seq2)):
        print(seq2[i], F[i])
    print()

    # back-tracking from max score
    output1 = []
    output2 = []
    row = m_row
    col = m_col
    while V[row][col] != 0:
        if V[row][col] == V[row-1][col-1]+align(seq1[col],seq2[row],letters,score_matrix):
            output1.insert(0, seq1[col])
            output2.insert(0, seq2[row])
            col -= 1
            row -= 1
        elif V[row][col] == F[row][col]:
            while F[row][col] == F[row-1][col]-s:
                output1.insert(0, "_")
                output2.insert(0, seq2[row])
                row -= 1
            output1.insert(0, "_")
            output2.insert(0, seq2[row])
            row -= 1
        elif V[row][col] == E[row][col]:
            while E[row][col] == E[row][col-1]-s:
                output1.insert(0, seq1[col])
                output2.insert(0, "_")
                col -= 1
            output1.insert(0, seq1[col])
            output2.insert(0, "_")
            col -= 1

    # print to console
    print("score = ", maximum)
    print("entries = ", entries)
    print()
    print(seq_names[0], end="")
    print(''.join([str(elem) for elem in output1]))
    print()
    print(seq_names[1], end="")
    print(''.join([str(elem) for elem in output2]))

    # output to file:
    f = open(outputfile_name, 'w')
    print("score = ", maximum, file=f)
    print("entries = ", entries, file=f)
    print("", file=f)
    print(seq_names[0], end="", file=f)
    print(''.join([str(elem) for elem in output1]), file=f)
    print("", file=f)
    print(seq_names[1], end="", file=f)
    print(''.join([str(elem) for elem in output2]), file=f)


def align(a, b,letters, matrix):
    row = letters.index(b)
    col = letters.index(a)
    value = matrix[row][col]
    return value


def readfile():
    # readfile format: parameter.txt input.txt output.txt
    filenames = input().split()
    paramfile_name = filenames[0]
    inputfile_name = filenames[1]
    outputfile_name = filenames[2]

    f = open(paramfile_name, "r")
    thresholdX = int(f.readline().split(";")[0])
    bandwidth = int(f.readline().split(";")[0])
    gap = int(f.readline().split(";")[0])
    indel = int(f.readline().split(";")[0])
    f.readline()
    f.readline()
    letters = f.readline().strip().split(" ")

    matrix = []
    f.readline()
    f.readline()
    lines = f.readlines()
    for line in lines:
        temp = []
        row = line.strip().split(" ")
        for x in row:
            if x != '':
                temp.append(int(x))
        if temp:
            matrix.append(temp)

    f = open(inputfile_name, "r")
    sequences = ""
    seq_names = []
    has_seq = False
    lines = f.readlines()
    for line in lines:
        if line[0] == '>':
            has_seq = True
            seq_names.append(line)
            continue
        if line == "\n" or line == "":
            has_seq = False
            sequences += "_"
        if has_seq:
            line = line.strip()
            sequences += line
    temp = sequences.split("_")
    seq1 = "_" + temp[0]
    seq2 = "_" + temp[1]
    for item in seq_names:
        item.strip()
    f.close()
    return matrix, letters, thresholdX, bandwidth, gap, indel, seq1, seq2, seq_names, outputfile_name


if __name__ == '__main__':
    main()