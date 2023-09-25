#该代码是根据GCG Wisconsin Package形式密码子表计算出给定序列的CAI值
'''
输入的密码子表为GCG Wisconsin Package形式密码子表，下载的网址为http://www.kazusa.or.jp/codon/

CAI计算原理
CAI = RSCUnn / RSCUmax

RSCU = RSCUFZ(某个密码子的个数)*n(该氨基酸密码子的种类数)/RSCUFM(该氨基酸总的密码子数)

具体可参考论文https://www.ncbi.nlm.nih.gov/pmc/articles/PMC340524/

'''
import re
#输入GCG Wisconsin Package形式密码子表的文件地址
infile = "ecoli.txt"
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
#计算的RSCU分母，即单个氨基酸的总密码子个数
RSCUFM = dict(zip(MM, values1))
#最高数量密码子的个数，只取占比最高的密码子，作为计算RSCUmax的分子
RSCUFZ = dict(zip(MM, values1))
#全部密码子的个数，作为计算RSCU的分子
RSCUFZ1 = dict(zip(encodings1, encodings2))
#记录RSCU分子的占比最高的密码子
RSCUFZA = dict(zip(MM, values2))
#记录计算的最高数量密码子的RSCU的值
RSCUNN = dict(zip(MM, values1))
#记录全部密码子的RSCU值
RSCUNN1 = dict(zip(encodings1, values3))
#记录氨基酸及氨基酸对应的密码子的个数
RSCUFZ_LIST = []
for i in range(42):
    RSCUFZ_LIST.append([])
for i in range(0,192):
    k = 0
    for j in MM:
        if encodings[i] == j:
            str1 = str(encodings[i+1])
            m = int(float(encodings[i+2]))
            RSCUFM[j] = RSCUFM[j] + m
            RSCUFZ_LIST[k].append(str1)
            RSCUFZ_LIST[k+1].append(m)
        k += 2
#记录计算RSCU里面的n值,即每个氨基酸对应的n值
mmm = [4, 2, 2, 4, 4, 6, 6, 2, 2, 1, 3, 4, 1, 3, 2, 2, 6, 2, 2, 2, 4]
#计算最高数量密码子的RSCU值
j = 0
k = 0
for i in MM:
    n = RSCUFZ_LIST[j+1].index(max(RSCUFZ_LIST[j+1]))
    RSCUFZA[i] = str(RSCUFZ_LIST[j][n])
    RSCUFZ[i] = int(RSCUFZ_LIST[j+1][n]) * mmm[k]
    j = j + 2
    k = k + 1
for i in MM:
    RSCUNN[i] = int(RSCUFZ[i]) / int(RSCUFM[i])
#记录密码子
AA =  'GATC'
triN = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]
#储存要除的分母,即密码子对应的氨基酸在列表MM的位置
mmm1 = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4,
        5, 5, 6, 6, 7, 7, 8, 8, 9, 10, 10, 10, 11, 11, 11, 11,
        12, 13, 14, 14, 13, 13, 15, 15, 16, 16, 17, 17, 6, 6, 6, 6,
        5, 5, 5, 5, 18, 18, 19, 19, 16, 16, 16, 16, 20, 20, 20, 20]
#存储要乘的RSCU的n值，即同义密码子个数
mmm2= [4, 4, 4, 4, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4,
       6, 6, 6, 6, 2, 2, 2, 2, 1, 3, 3, 3, 4, 4, 4, 4,
       1, 3, 2, 2, 3, 3, 2, 2, 6, 6, 2, 2, 6, 6, 6, 6,
       6, 6, 6, 6, 2, 2, 2, 2, 6, 6, 6, 6, 4, 4, 4, 4]
#计算全部密码子的RSCU值
for i in range(64):
    RSCUNN1[triN[i]] = (RSCUFZ1[triN[i]]/RSCUFM[MM[mmm1[i]]]) * int(mmm2[i])
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
condonAA = dict(zip(encodings1, encodings0))
def CAIC(sequence):
    #记录CAI值
    CAI0 = 1
    #记录序列密码子个数
    n = 0
    for i in range(0, len(sequence), 3):
        if i+3 < len(sequence):
            sequencesTTT = sequence[i : i + 3].upper()
            CAI0 = CAI0 *(RSCUNN1[sequencesTTT] / RSCUNN[condonAA[sequencesTTT]])
            n += 1
    CAIN = CAI0 ** (1/n)
    return CAIN
CAIcounter = []
for sequence in sequences:
    CAIcounter.append(CAIC(sequence))
for i in range(len(f1)):
    print('序列：', f1[i][0], '在物种文件',infile,'的CAI=', CAIcounter[i])
