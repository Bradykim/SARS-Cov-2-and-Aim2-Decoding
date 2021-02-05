from compsci260lib import *   # this will also import the sys and re modules


def run_orfs():
    """"
    Report the number of ORFs if the minimum ORF length is 10 amino acids, 
    40 amino acids, or 60 amino acids.  Also report the average length (in 
    amino acids) of the identified ORFs for the three cases.
    """

    sars_cov2_fasta_path = 'sars_cov2_wu.fasta'
    sars_cov2_temp = get_fasta_dict(sars_cov2_fasta_path)
    sars_cov2_genomic = sars_cov2_temp['MN908947.3'].upper()
    finalOrfs10 = find_orfs(sars_cov2_genomic, 10)
    numOrfs10 = summarize_orfs(finalOrfs10)[0]
    avglengthOrfs10 = summarize_orfs(finalOrfs10)[1]
    finalOrfs40 = find_orfs(sars_cov2_genomic, 40)
    numOrfs40 = summarize_orfs(finalOrfs40)[0]
    avglengthOrfs40 = summarize_orfs(finalOrfs40)[1]
    finalOrfs60 = find_orfs(sars_cov2_genomic, 60)
    numOrfs60 = summarize_orfs(finalOrfs60)[0]
    avglengthOrfs60 = summarize_orfs(finalOrfs60)[1]
    report = f'Orf length of 10 amino acids: number of orfs is {numOrfs10} and the average length of the orfs is {avglengthOrfs10}. \n' \
             f'Orf length of 40 amino acids: number of orfs is {numOrfs40} and the average length of the orfs is {avglengthOrfs40}. \n' \
             f'Orf length of 60 amino acids: number of orfs is {numOrfs60} and the average length of the orfs is {avglengthOrfs60}.'
    print(report)
    print(find_orfs(sars_cov2_genomic,0))
    return report


def summarize_orfs(orfs):
    """Summarize ORFs identified from the find_orfs procedure as a count of the
    number of found orfs and the average length of the found ORFs (in amino
    acids)

    Args:
        orfs (list): a list of dictionaries of found ORFs

    Returns:
        tuple: (The number of ORFs found (int), Average ORF length (float))
    """

    numORFs = len(orfs)
    total =0

    for item in orfs:
        total+=item['aalength']

    averageORFs = total/numORFs
    return (numORFs, averageORFs)

def orfend(seq, start):
    """
    This is a function for when the start codon AUG is found, to loop through the
    sequence until a stop codon is found while gathering all necessary information. 
    
    The function takes parameters: a genomic sequence, the index of the G of 'AUG'
    """
    finalIndex = start
    finalseq = ''
    for i in range(start, len(seq), 3):
        codon = seq[i: i + 3]
        finalIndex += 3
        if (codon == 'UAG') or (codon == 'UGA') or (codon == 'UAA'):
            finalseq = codon
            break

    tuple = (finalseq, finalIndex)
    return tuple

def internalOrfs(orflist):
    for i in range(len(orflist)):
        for j in range(i + 1, len(orflist)):
            if orflist[i]['stop'] == orflist[j]['stop']:
                if orflist[i]['start'] < orflist[j]['start']:
                    orflist[j]['remove'] = 1
                else:
                    orflist[i]['remove'] = 1
    orflist[:] = [x for x in orflist if not 'remove' in x.keys()]
    return orflist


def find_orfs(seq, min_length_aa):
    """This is a function for finding sufficiently long ORFs in all reading
    frames in a sequence of DNA or RNA.  By default, the sequence is assumed
    to be single-stranded.

    The function takes as input parameters: a string representing a genomic
    sequence, the minimum length (in amino acids) for an ORF before it will be
    returned (which defaults to 0).

    Args:
        seq (str): a genomic sequence
        min_length_aa (int): minimum length of found ORFs in amino acids

    Returns:
        list: of dictionaries with information on each ORF found.

    Where each ORF found is represented by a dictionary with
    the following keys:
        frame (int): The nucleotide offset in which the ORF was found. (Must be
        0, 1, or 2)
        stop (int): the nucleotide position of the end of the ORF
        start (int): the nucleotide position of the start of the ORF
        stopcodon (str): the nucleotide triplet of the stop codon
        nlength (int): the length (in nucleotides) of the ORF
        strand (str): the strand of the found ORF (Must be '+' or '-')

    A valid return list may look something like this:
    [
        {
            'frame': 0,
            'stop': 13413,
            'aalength': 4382,
            'start': 265,
            'stopcodon': 'UAA',
            'nlength': 13149,
            'strand': '+'
        },
        {
            'frame': 0,
            'stop': 27063,
            'aalength': 221,
            'start': 26398,
            'stopcodon': 'UAA',
            'nlength': 666,
            'strand': '-'
        }
    ]
    """
    orfs = []
    AUGIiterator = re.finditer('AUG', seq)
    for match in AUGIiterator:
        finalseq = orfend(seq, match.start())
        length = finalseq[1] - match.start()
        aalength = int(length/3)
        if(aalength>= min_length_aa) and not (finalseq[0] == ''):
            orfs.append(
                {'frame': match.start() % 3, 'stop': finalseq[1], 'aalength': aalength, 'start': match.start(),
                 'stopcodon': finalseq[0], 'nlength': length, 'strand': '+'})
    orfsfinal = internalOrfs(orfs)
    return orfsfinal


if __name__ == '__main__':
    """Run run_orfs(). Do not modify this code"""
    run_orfs()
