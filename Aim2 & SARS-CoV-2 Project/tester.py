from compsci260lib import *
# def location(startaim2, endaim2, startrestriction, endrestriction ):
#     if (startrestriction>=startaim2) and (endrestriction<=endaim2):
#         return 'within'
#     elif endaim2 <= startrestriction:
#         return 'downstream'
#     else:
#         return 'upstream'
#
# aim2_beg = 1000
# aim2_end = 1740
# aim2_fasta_path = 'aim2_plus_minus_1kb.fasta'
# aim2_genomic = get_fasta_dict(aim2_fasta_path)
#
#
# # r_enzymes = {
# #     'BamHI': 'GGATCC',
# #     'BstYI': r'[A|G]GATC[C|T]',
# #     'SalI': 'GTCGAC',
# #     'SpeI': 'ACTAGT',
# #     'SphI': 'GCATGC',
# #     'StyI': 'CCATGG'
# # }
# #
# #
# # r_enzyme_sites = {"BamHI": [], "BstYI": [], "SpeI": [],
# #                       "SphI": [], "StyI": [], "SalI": []}
#
# # genomicSequence = aim2_genomic['aim2'].upper()
# # BamHIiterator = re.finditer(r_enzymes['BamHI'], genomicSequence)
# # for match in BamHIiterator:
# #     currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group(), 'location': location(aim2_beg,aim2_end,match.start(),match.end())}
# #     r_enzyme_sites['BamHI'].append(currentDict)
# # BstYIiterator = re.finditer(r_enzymes['BstYI'], genomicSequence)
# # for match in BstYIiterator:
# #     currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group(), 'location': location(aim2_beg,aim2_end,match.start(),match.end())}
# #     r_enzyme_sites['BstYI'].append(currentDict)
# # SpeIiterator = re.finditer(r_enzymes['SpeI'], genomicSequence)
# # for match in SpeIiterator:
# #     currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group(), 'location': location(aim2_beg,aim2_end,match.start(),match.end())}
# #     r_enzyme_sites['SpeI'].append(currentDict)
# # SphIiterator = re.finditer(r_enzymes['SphI'], genomicSequence)
# # for match in SphIiterator:
# #     currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group(), 'location': location(aim2_beg,aim2_end,match.start(),match.end())}
# #     r_enzyme_sites['SphI'].append(currentDict)
# # StyIiterator = re.finditer(r_enzymes['StyI'], genomicSequence)
# # for match in StyIiterator:
# #     currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group(), 'location': location(aim2_beg,aim2_end,match.start(),match.end())}
# #     r_enzyme_sites['StyI'].append(currentDict)
# # SalIiterator = re.finditer(r_enzymes['SalI'], genomicSequence)
# # for match in SalIiterator:
# #     currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group(), 'location': location(aim2_beg,aim2_end,match.start(),match.end())}
# #     r_enzyme_sites['SalI'].append(currentDict)
# # print(r_enzyme_sites)
# covid_fasta_path = 'sars_cov2_wu.fasta'
# covid_genomic = get_fasta_dict(covid_fasta_path)
# genomicSequence = covid_genomic['MN908947.3'].upper()
# def orfend(seq, start):
#     finalIndex = 0
#     finalseq = ''
#     for i in range(start, len(seq)):
#         codon = seq[i * 3: i * 3 + 3]
#         if (codon == 'UAG') or (codon == 'UGA') or (codon == 'UAA'):
#             finalseq = codon
#             break
#         finalIndex += 3
#     tuple = (finalseq,start+finalIndex)
#     return tuple
#
# def internalOrfs(orflist):
#     for i in range(len(orflist)):
#         for j in range(i+1,len(orflist)):
#             if orflist[i]['stop'] == orflist[j]['stop']:
#                 if(orflist[i]['start']<orflist[j]['start']):
#                     orflist[j]['remove'] =1
#                 else:
#                     orflist[i]['remove'] =1
#     orflist[:] = [x for x in orflist if not 'remove' in x.keys()]
#     return orflist
#
# def orfssss():
#     numORFs = len(orfs1())
#     total = 0
#
#     for item in orfs1():
#         total += item['aalength']
#
#     averageORFs = total / numORFs
#     return (numORFs, averageORFs)
#
# def orfs1():
#     orfs =[]
#     AUGIiterator = re.finditer('AUG', genomicSequence)
#     for match in AUGIiterator:
#         finalseq = orfend(genomicSequence, match.end())
#         length = finalseq[1]-match.start()
#         orfs.append({'frame': match.start()%3, 'stop': finalseq[1], 'aalength': int(length/3),'start': match.start(),
#                      'stopcodon': finalseq[0],'nlength':length, 'strand':'+'})
#     orfsfinal = internalOrfs(orfs)
#     return orfs;
#
#
# print(orfssss())
# print(orfs1())
# print(genomicSequence)
import math

seq='AAAAUGAAAUGA'
start = 3
finalIndex = start
finalseq = ''
for i in range(start, len(seq),3):
    codon = seq[i : i + 3]
    print(codon)
    finalIndex+=3
    if (codon == 'UAG') or (codon == 'UGA') or (codon == 'UAA'):
        finalseq = codon
        break


tuple = (finalseq, finalIndex)

print (tuple)

