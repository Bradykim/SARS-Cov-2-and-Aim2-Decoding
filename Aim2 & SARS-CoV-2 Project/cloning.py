    
from compsci260lib import *   # this will also import the sys and re modules


def run_cloning():
    """Function demonstrating the procedure of extracting a gene and inserting
    into a plasmid using restriction enzymes."""

    aim2_fasta_path = 'aim2_plus_minus_1kb.fasta'
    pRS304_fasta_path = 'pRS304.fasta'

    # Read in the Aim2 genomic sequence from the fasta file along with its
    # upstream and downstream regions.
    #
    # Your code goes here
    #
    aim2temp = get_fasta_dict(aim2_fasta_path)
    aim2_genomic = aim2temp['aim2'].upper()

    # Store the beginning and end locations of the Aim2 gene as Python indices.
    aim2_beg = 1001  # Your value goes here
    aim2_end = 1741  # Your value goes here
    aim2_length = (aim2_end-aim2_beg+1)
    aim2_aminolength = int(aim2_length/3)
    # Report start, end, length in nucleotides, and length in amino acids.
    print(f'The starting location of the gene is {aim2_beg}.')
    print(f'The ending position of the gene is {aim2_end}.')
    print(f'The length of the gene in nucleotides is {aim2_length}.')
    print(f'The length of the gene in amino acids is {aim2_aminolength}.')
    # Define regular expression terms for each restriction enzyme
    r_enzymes = get_restriction_enzymes_regex()

    # Store coordinates of restriction sites found upstream, downstream, and
    # within the aim2 gene
    r_enzyme_sites = find_aim2_restriction_enzymes(
        aim2_beg, aim2_end, aim2_genomic)

    # Report the found restriction enzyme sites
    #
    restrictionEnzymes = find_aim2_restriction_enzymes(aim2_beg,aim2_end,aim2_genomic)
    print(f'\nBamHI: No enzymes found.\n')

    BstYIlist = restrictionEnzymes['BstYI']
    print(f'BstYI: ')
    for item in BstYIlist:
        start = item.get('start')
        end = item.get('end')
        sequence = item.get('sequence')
        location = item.get('location')
        print(f'enzyme start is {start}, end is {end}, sequence is {sequence}, location is {location}')

    print(f'\nSpeI: No enzymes found.\n')

    SphIlist = restrictionEnzymes['SphI']
    print(f'SphI: ')
    for item in SphIlist:
        start = item.get('start')
        end = item.get('end')
        sequence = item.get('sequence')
        location = item.get('location')
        print(f'enzyme start is {start}, end is {end}, sequence is {sequence}, location is {location}')

    print(f'\nStyI: No enzymes found.\n')

    SalIlist = restrictionEnzymes['SalI']
    print(f'SalI: ')
    for item in SalIlist:
        start = item.get('start')
        end = item.get('end')
        sequence = item.get('sequence')
        location = item.get('location')
        print(f'enzyme start is {start}, end is {end}, sequence is {sequence}, location is {location} \n')

    #

    # Read in the pRS304 plasmid sequence and find restriction sites
    # in the entire plasmid
    prs304temp = get_fasta_dict(pRS304_fasta_path)
    prs304_genomic = prs304temp['pRS304'].upper()
    all_prs304_sites = find_pRS304_restriction_sites(prs304_genomic)

    # Report relevant summary information
    mcs_start =  1886
    mcs_end = 1993
    mcs_enzyme_sites = report_pRS304_MCS_sites(all_prs304_sites, mcs_start,
        mcs_end)
    bamout = mcs_enzyme_sites['BamHI']['outside']
    bamin = mcs_enzyme_sites['BamHI']['inside']
    print(f'BamHI: outside is {bamout}, inside is {bamin}')

    bstout = mcs_enzyme_sites['BstYI']['outside']
    bstin = mcs_enzyme_sites['BstYI']['inside']
    print(f'BstYI: outside is {bstout}, inside is {bstin}')

    speout = mcs_enzyme_sites['SpeI']['outside']
    spein = mcs_enzyme_sites['SpeI']['inside']
    print(f'SpeI: outside is {speout}, inside is {spein}')

    sphout = mcs_enzyme_sites['SphI']['outside']
    sphin = mcs_enzyme_sites['SphI']['inside']
    print(f'SphI: outside is {sphout}, inside is {sphin}')

    styout = mcs_enzyme_sites['StyI']['outside']
    styin = mcs_enzyme_sites['StyI']['inside']
    print(f'StyI: outside is {styout}, inside is {styin}')

    salout = mcs_enzyme_sites['SalI']['outside']
    salin = mcs_enzyme_sites['SalI']['inside']
    print(f'SalI: outside is {salout}, inside is {salin}')
    # Extract aim2 gene and insert into the plasmid, report its length
    #
    # Your code goes here
    #


def get_restriction_enzymes_regex():
    """Returns a dictionary of restriction enzyme regular expressions for
    searching in genomic sequences.

    This function should be used for find_aim2_restriction_enzymes and
    find_pRS304_MCS_restriction_sites.
    """

    r_enzymes = {
        "BamHI": "GGATCC",
        "BstYI": r'[A|G]GATC[C|T]',
        "SalI": "GTCGAC",
        "SpeI": "ACTAGT",
        "SphI": "GCATGC",
        "StyI": "CC[A|T][A|T]GG"
    }
    return r_enzymes

def location(startaim2, endaim2, startrestriction, endrestriction ):
    """
    Returns the location of the restriction enzyme relative to the
    location of the Aim2 gene.

    Args:
        startaim2 (int)
        endaim2 (int)
        startrestriction (int)
        endrestriction (int)

    Will return values of:
    'within': if enzyme is located in Aim2
    'upstream': if enzyme is located before Aim2 in the gene
    'downstream': if enzyme is located after Aim2 in the gene

    """
    if (startrestriction>=startaim2) and (endrestriction<=endaim2):
        return 'within'
    elif endaim2 <= startrestriction:
        return 'downstream'
    else:
        return 'upstream'



def find_aim2_restriction_enzymes(aim2_beg, aim2_end, aim2_genomic):

    """Finds the restriction enzyme sites in the aim2 genomic sequence. Stored
    as a dictionary of lists. Each restriction enzyme key corresponds to a list
    of dictionaries containing relevant information for a found site.

    Args:
        aim2_beg (int)
        aim2_end (int)
        aim2_genomic (str): genomic sequence to search for sites in

    Each found site will have defined keys:
        sequence (str): the nucleotide sequence matched in the genome
        start (int): the start nucleotide position in the genome
        end (int): the ending position
        location (str): the position of the site relative the aim2 gene.
            (must be "upstream", "downstream", or "within")

    A valid returned dictionary may look like this:
    {
        'BamHI' : [
            {
                'start': 10,
                'end': 15,
                'sequence': 'GGATCC',
                'location': 'upstream'
            },
            {
                'start': 100,
                'end': 105,
                'sequence': 'GGATCC',
                'location': 'downstream'
            }
        ],
        'BstYI' : [
            {
                'start': 30,
                'end': 35,
                'sequence': 'AGATCC',
                'location': 'within'
            }
        ]
    }
    """

    if type(aim2_genomic) is not str:
        raise TypeError(f"the argument 'aim2_genomic' must be a string. "
                        "(received {type(aim2_genomic)})")

    # load restriction enzyme regular expressions
    r_enzymes = get_restriction_enzymes_regex()

    # Create dictionary to store coordinates
    r_enzyme_sites = {"BamHI": [], "BstYI": [], "SpeI": [],
                      "SphI": [], "StyI": [], "SalI": []}

    #
    # Your code goes here
    #
    r_enzymes = get_restriction_enzymes_regex()
    genomicSequence = aim2_genomic

    BamHIiterator = re.finditer(r_enzymes['BamHI'], genomicSequence)
    for match in BamHIiterator:
        if isPalinrome(match.group()):
            currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group(), 'location': location(aim2_beg,aim2_end,match.start(),match.end())}
            r_enzyme_sites['BamHI'].append(currentDict)

    BstYIiterator = re.finditer(r_enzymes['BstYI'], genomicSequence)
    for match in BstYIiterator:
        if isPalinrome(match.group()):
            currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group(), 'location': location(aim2_beg,aim2_end,match.start(),match.end())}
            r_enzyme_sites['BstYI'].append(currentDict)

    SpeIiterator = re.finditer(r_enzymes['SpeI'], genomicSequence)
    for match in SpeIiterator:
        if isPalinrome(match.group()):
            currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group(), 'location': location(aim2_beg,aim2_end,match.start(),match.end())}
            r_enzyme_sites['SpeI'].append(currentDict)

    SphIiterator = re.finditer(r_enzymes['SphI'], genomicSequence)
    for match in SphIiterator:
        if isPalinrome(match.group()):
            currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group(), 'location': location(aim2_beg,aim2_end,match.start(),match.end())}
            r_enzyme_sites['SphI'].append(currentDict)

    StyIiterator = re.finditer(r_enzymes['StyI'], genomicSequence)
    for match in StyIiterator:
        if isPalinrome(match.group()):
            currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group(), 'location': location(aim2_beg,aim2_end,match.start(),match.end())}
            r_enzyme_sites['StyI'].append(currentDict)

    SalIiterator = re.finditer(r_enzymes['SalI'], genomicSequence)
    for match in SalIiterator:
        if isPalinrome(match.group()):
            currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group(), 'location': location(aim2_beg,aim2_end,match.start(),match.end())}
            r_enzyme_sites['SalI'].append(currentDict)
    return r_enzyme_sites


def isPalinrome(str):
    if(str==reverse_complement(str)):
        return True
    else:
        return False
def report(list, mcs_start, mcs_end):
    outside =0
    inside =0
    iter=0
    for i in list:
        if location(list[iter]['start'],list[iter]['end'],mcs_start,mcs_end) == 'within':
            inside+=1
        else:
            outside+=1
        iter+=1
    tuple =(outside, inside)
    return tuple
def report_pRS304_MCS_sites(p_enzyme_sites, mcs_start, mcs_end):
    """
    For each restriction enzyme,

    Report how often that enzyme cuts the plasmid outside the MCS, how often
    it cuts the plasmid inside the MCS, and relevant details about any sites
    located inside the MCS.
    """
    r_enzymes = {"BamHI": {'outside': 0, 'inside': 0}, "BstYI": {'outside': 0, 'inside': 0}, "SpeI": {'outside': 0, 'inside': 0},
                      "SphI": {'outside': 0, 'inside': 0}, "StyI": {'outside': 0, 'inside': 0}, "SalI": {'outside': 0, 'inside': 0}}

    restrictenzymes = list(p_enzyme_sites.keys())
    for x in restrictenzymes:
        if x == 'BamHI':
            tup = report(p_enzyme_sites['BamHI'], mcs_start, mcs_end)
            r_enzymes['BamHI']['outside'] = tup[0]
            r_enzymes['BamHI']['inside'] = tup[1]
        if x == 'BstYI':
            tup = report(p_enzyme_sites['BstYI'], mcs_start, mcs_end)
            r_enzymes['BstYI']['outside'] = tup[0]
            r_enzymes['BstYI']['inside'] = tup[1]
        if x == 'SpeI':
            tup = report(p_enzyme_sites['SpeI'], mcs_start, mcs_end)
            r_enzymes['SpeI']['outside'] = tup[0]
            r_enzymes['SpeI']['inside'] = tup[1]
        if x == 'SphI':
            tup = report(p_enzyme_sites['SphI'], mcs_start, mcs_end)
            r_enzymes['SphI']['outside'] = tup[0]
            r_enzymes['SphI']['inside'] = tup[1]
        if x == 'StyI':
            tup = report(p_enzyme_sites['StyI'], mcs_start, mcs_end)
            r_enzymes['StyI']['outside'] = tup[0]
            r_enzymes['StyI']['inside'] = tup[1]

        if x == 'SalI':
            tup = report(p_enzyme_sites['SalI'], mcs_start, mcs_end)
            r_enzymes['SalI']['outside'] = tup[0]
            r_enzymes['SalI']['inside'] = tup[1]



    return  r_enzymes


def find_pRS304_restriction_sites(prs304_genomic):
    """Finds the restriction sites in pRS304 genomic sequence. Stored as a
    dictionary of lists. Each restriction enzyme key corresponds to a list of
    dictionaries containing relevant information for a found site.


    This code will be similar (and simpler) code
    found in `find_aim2_restriction_enzymes`. So it may be helpful to use code
    you wrote there.

    Args:
        prs304_genomic (str): genomic sequence to search for sites in

    Each found site will have defined keys:
        sequence (str): the nucleotide sequence matched
        start (int): the start nucleotide position in the plasmid
        end (int): the ending position in the plasmid

    A valid returned dictionary may look like this:
    {
        'BamHI' : [
            {
                'start': 10,
                'end': 15,
                'sequence': 'GGATCC'
            },
            {
                'start': 100,
                'end': 105,
                'sequence': 'GGATCC'
            }
        ],
        'BstYI' : [
            {
                'start': 30,
                'end': 35,
                'sequence': 'AGATCC'
            }
        ]
    }
    """

    if type(prs304_genomic) is not str:
        raise TypeError(f"the argument 'prs304_genomic' must be a string. "
                        "(received type {type(prs304_genomic)})")

    # load restriction enzyme regular expressions
    r_enzymes = get_restriction_enzymes_regex()

    # Search for restriction sites in the pRS304 MCS.
    p_enzyme_sites = {"BamHI": [], "BstYI": [], "SpeI": [],
                      "SphI": [], "StyI": [], "SalI": []}

    genomicSequence = prs304_genomic

    BamHIiterator = re.finditer(r_enzymes['BamHI'], genomicSequence)
    for match in BamHIiterator:
        currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group()}
        p_enzyme_sites['BamHI'].append(currentDict)
    BstYIiterator = re.finditer(r_enzymes['BstYI'], genomicSequence)
    for match in BstYIiterator:
        currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group()}
        p_enzyme_sites['BstYI'].append(currentDict)
    SpeIiterator = re.finditer(r_enzymes['SpeI'], genomicSequence)
    for match in SpeIiterator:
        currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group()}
        p_enzyme_sites['SpeI'].append(currentDict)
    SphIiterator = re.finditer(r_enzymes['SphI'], genomicSequence)
    for match in SphIiterator:
        currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group()}
        p_enzyme_sites['SphI'].append(currentDict)
    StyIiterator = re.finditer(r_enzymes['StyI'], genomicSequence)
    for match in StyIiterator:
        currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group()}
        p_enzyme_sites['StyI'].append(currentDict)
    SalIiterator = re.finditer(r_enzymes['SalI'], genomicSequence)
    for match in SalIiterator:
        currentDict = {'start': match.start(), 'end': match.end(), 'sequence': match.group()}
        p_enzyme_sites['SalI'].append(currentDict)

    return p_enzyme_sites


if __name__ == '__main__':
    """Run run_cloning(). Do not modify this code"""
    run_cloning()
