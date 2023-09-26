#该代码是根据GCG Wisconsin Package形式密码子表计算出给定序列的CAI值
'''
输入的密码子表为GCG Wisconsin Package形式密码子表，下载的网址为http://www.kazusa.or.jp/codon/

CAI计算原理
CAI = RSCUnn / RSCUmax
由于
RSCU = RSCUFZ(某个密码子的个数)*n(该氨基酸密码子的种类数)/RSCUFM(该氨基酸总的密码子数)
具体可参考论文https://www.ncbi.nlm.nih.gov/pmc/articles/PMC340524/
则可进一步简化为
CAI = RSCUnn / RSCUmax = Nij / Nmax ** (1/n)
即该位置密码子的该密码子在该物种的数量除以该位置该密码子编码的氨基酸在该物种对应的最多的密码子的数量再开方，开该序列密码子数量的次方

'''
import re
#输入GCG Wisconsin Package形式密码子表的文件地址
infile = "test.txt"
#输入序列文件的文件地址
infile1 = "testseq.txt"
#读取密码子表
with open(infile) as f:
    records = f.read()
records = records.split()
encodings = []
encodings0 = []
encodings1 = []
encodings2 = []
for i in range(0,320,5):
    encodings.append(records[i])
    encodings.append(records[i+1])
    encodings.append(records[i+2])
    encodings0.append(records[i])
    encodings1.append(records[i+1])
    encodings2.append(int(float(records[i+2])))
#记录氨基酸
MM = ['Gly', 'Glu', 'Asp', 'Val', 'Ala','Arg', 'Ser', 'Lys', 'Asn', 'Met',
    'Ile', 'Thr', 'Trp', 'End', 'Cys','Tyr', 'Leu', 'Phe', 'Gln', 'His', 'Pro']
values1 = [0] * 21
values2 = ['QQQ'] * 21
values3 = [0] * 64

#最高数量密码子的个数，只取占比最高的密码子，作为计算CAI的分母
CAINmax = dict(zip(MM, values1))
#全部密码子每个密码子的数量，作为计算CAI的分子
CAINij = dict(zip(encodings1, encodings2))
#记录每个密码子的CAI值
CAINN = dict(zip(encodings1, values3))
#记录氨基酸对应占比最高的密码子
CAINmaxcodon = dict(zip(MM, values2))
#记录氨基酸及氨基酸对应的密码子的个数，作为找最高密码子的缓冲
CAINmax_LIST = []
for i in range(42):
    CAINmax_LIST.append([])
for i in range(0,192):
    k = 0
    for j in MM:
        if encodings[i] == j:
            str1 = str(encodings[i+1])
            m = int(float(encodings[i+2]))
            CAINmax_LIST[k].append(str1)
            CAINmax_LIST[k+1].append(m)
        k += 2
#记录每个氨基酸最高数量密码子的密码子
j = 0
k = 0
for i in MM:
    n = CAINmax_LIST[j+1].index(max(CAINmax_LIST[j+1]))
    CAINmaxcodon[i] = str(CAINmax_LIST[j][n])
    CAINmax[i] = CAINmax_LIST[j+1][n]
    j = j + 2
    k = k + 1
#记录密码子
AA =  'GATC'
triN = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]
condonAA = dict(zip(encodings1, encodings0))
#计算全部密码子的CAI值
for i in range(64):
    CAINN[triN[i]] = CAINij[triN[i]]/CAINmax[condonAA[triN[i]]]
#读取序列
def get_fasta(file):
    text = []
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    for fasta in records:
        array = fasta.split('\n')
        header, sequence = array[0].split()[0], re.sub('[^ACGTU-]', '-', ''.join(array[1:]).upper())
        sequence = re.sub('U', 'T', sequence)
        text.append([header, sequence])
    return text

f1 = get_fasta(infile1)
sequences = []
for i in range(len(f1)):
    sequences.append(f1[i][1])
#密码子转换为氨基酸的表
def CAIC(sequence):
    #记录CAI值
    CAI0 = 1
    #记录序列密码子个数
    n = 0
    for i in range(0, len(sequence), 3):
        if i+3 < len(sequence):
            sequencesTTT = sequence[i : i + 3].upper()
            CAI0 = CAI0 * CAINN[sequencesTTT]
            n += 1
    CAIN = CAI0 ** (1/n)
    return CAIN
CAIcounter = []
for sequence in sequences:
    CAIcounter.append(CAIC(sequence))
for i in range(len(f1)):
    print('序列：', f1[i][0], '在密码子表文件',infile,'的CAI=', CAIcounter[i])


