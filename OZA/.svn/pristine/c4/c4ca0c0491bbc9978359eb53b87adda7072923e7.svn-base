import numpy as np

# Lezen van Deltares TEKAL blokken en tabel-interpolatie
def getTEKALLine(fn):
    while True:
        line = fn.readline()
        if not(line):
             return(line)
        if (line[0]!="*"):
             return(line)

def readTEKAL(fnin,blklbl,**kwargs):
# Read selected TEKAL block in a list of numpy 1d-arrays
    try:
        each = int(kwargs['skip'])+1
    except:
        each = 1
    with open(fnin,"r") as fn:
        cnt = 0
        while True:
            line = getTEKALLine(fn)
            if not(line):
                raise (Exception("Label "+blklbl+" not found!"))
            if line.strip()==blklbl:
                line = getTEKALLine(fn)
                try:
                    nrowcol = line.split()
                    nrow, ncol = int(nrowcol[0]), int(nrowcol[1])
                    break
                except:
                    raise (Exception("Block "+blklbl+" corrupt!"))
        data =[]   # list of numpy arrays
        for icol in range(ncol):
            data.append(np.array([]))
        for irow in range(nrow):
            line = getTEKALLine(fn)
            if not(line):
                break
            strings = line.split()
            if cnt%each == 0:
                for icol in range(ncol):
                    data[icol] = np.append(data[icol],float(strings[icol]))
            cnt += 1
    return data


def listTEKALBlocks(fnin):
    blocks = {}
    with open(fnin,"r") as fn:
        pos = 0
        while True:
            line = fn.readline()
            pos += 1
            if not line:
                break
            if line[0] in '*#':
                continue
            lbl = line.strip()  
            blocks[lbl] = {}
            line = fn.readline()
            pos += 1
            while line[0] in '*#': 
                line = fn.readline()
                pos += 1
            nrow, ncol = [int(substr.strip()) for substr in line.split()]
            blocks[lbl]['nrow'] = nrow 
            blocks[lbl]['ncol'] = ncol
            blocks[lbl]['linenr'] = pos     # first line following the row column specification
            for i in range(nrow):
                line = fn.readline()        # skip nrow lines
                pos += 1
    return blocks

