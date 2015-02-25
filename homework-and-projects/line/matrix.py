def new_matrix(rows = 4, cols = 4):
    m = []
    for c in range(cols):
        m.append([])
        for r in range(rows):
            m[c].append(0)
    return m

def print_matrix(matrix, in_pairs = 0):
    text = ""
    if(not in_pairs):
        for c in matrix:
            for r in c:
                text += str(r) + " "
            text += "\t"
    else:
        for index, item in enumerate(matrix):
            text += str(item) + "\t"
            if(index%2):
                text += "\n"        
    print text
    return
