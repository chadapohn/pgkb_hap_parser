import pandas as pd
import numpy as np
from pandas.testing import assert_frame_equal

from pgkb_hap_parser.parser import parse_alleles

def test_parse_alleles():
    """
    Output format:
            name    alleles
    0  Reference      G,C
    1   c.520C>T      A,C
    2  c.3257G>A      G,T
    """

    #################################################################
    # 1. Test parsing simple allele definition table 
    #################################################################
    obs = parse_alleles(
        'CACNA1S',
        pd.DataFrame(
            data=[
                ['Reference', 'G', 'C'],
                ['c.520C>T', 'A', np.NaN], 
                ['c.3257G>A', np.NaN, 'T']
            ]
        ),
        2
    )
    exp = pd.DataFrame(
        data={
            'name': ['Reference', 'c.520C>T', 'c.3257G>A'],
            'alleles': ['G,C', 'A,C', 'G,T']
        }
    )
    assert_frame_equal(obs, exp)
    
    #################################################################
    # 2. Test parsing allele definition table with notes
    #################################################################
    obs = parse_alleles(
        'CFTR',
        pd.DataFrame(
            data=[
                ['Reference', 'G', 'C', 'C', 'G', 'C'],
                ['F508del', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN],
                ['2789+5G->A', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN],
                ['3272-26A->G', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN], 
                ['3849+10kbC->T', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN],
                ['711+1G->T', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN], 
                ['A1067T', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN], 
                ['A455E', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN], 
                ['D110E', np.NaN, np.NaN, np.NaN, np.NaN, 'A'], 
                ['D110H', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN], 
                ['D1152H', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN], 
                ['D1270N', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN], 
                ['D579G', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN], 
                ['E193K', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN], 
                ['E56K', 'A', np.NaN, np.NaN, np.NaN, np.NaN],  
                [np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN],
                ['NOTES:', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN],
                ['rsID=RefSNP accession ID number (http://www.ncbi.nlm.nih.gov/snp/)', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN],
                ['All variants are annotated to the positive chromosomal strand.', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN],
                ['Class is defined in PMID: 20932301.', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN],
                ['The F508del CFTR variant can result from a CTT deletion at cDNA position NM_000492.3.c.1521_1523 (rs113993960) or a TCT deletion at cDNA position NM_000492.3.c.1520_1522 (rs199826652). Both result in the same sequence change: ATC ATC TTT GGT GTT > ATC ATT GGT GTT, corresponding to a deletion of Phe at amino acid position 508. Here we include rs113993960 (deletion CTT), the cDNA reference position name that is referred to on the CFTR1 website (c.1521_1523delCTT). Rs199826652 (deletion TCT)  is more likely to be called in sequencing data due to the left justification of indels; hence this has a minor allele frequency from 1000 genomes.', np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]
            ]
        ),
        5
    )
    exp = pd.DataFrame(
        data={
            'name': ['Reference', 'F508del', '2789+5G->A', '3272-26A->G', '3849+10kbC->T', '711+1G->T', 'A1067T', 'A455E', 'D110E', 'D110H', 'D1152H', 'D1270N', 'D579G', 'E193K', 'E56K'],
            'alleles': ['G,C,C,G,C', 'G,C,C,G,C', 'G,C,C,G,C', 'G,C,C,G,C', 'G,C,C,G,C', 'G,C,C,G,C', 'G,C,C,G,C', 'G,C,C,G,C', 'G,C,C,G,A', 'G,C,C,G,C', 'G,C,C,G,C', 'G,C,C,G,C', 'G,C,C,G,C', 'G,C,C,G,C', 'A,C,C,G,C']
        }
    )
    assert_frame_equal(obs, exp)
    
    #################################################################
    # 3. Test parsing allele definition table with functions
    #################################################################
    obs = parse_alleles(
        'G6PD',
        pd.DataFrame(
            data=[
                ['B (wildtype)', 'IV/Normal', 'G', 'C', 'C', 'T', 'T'],
                ['No name', 'Not reported/Unknown', 'A', np.NaN, np.NaN, np.NaN, np.NaN], 
                ['Sinnai', 'III/Deficient', np.NaN, 'A', np.NaN, np.NaN, np.NaN],
                ['Lages', 'III/Deficient', np.NaN, np.NaN, 'T', np.NaN, np.NaN],
                ['Gaohe', 'III/Deficient', np.NaN, np.NaN, np.NaN, 'C', np.NaN]
            ]
        ),
        5
    )
    exp = pd.DataFrame(
        data={
            'name': ['B (wildtype)', 'No name', 'Sinnai', 'Lages', 'Gaohe'],
            'alleles': ['G,C,C,T,T', 'A,C,C,T,T', 'G,A,C,T,T', 'G,C,T,T,T', 'G,C,C,C,T']
        }
    )
    assert_frame_equal(obs, exp)
