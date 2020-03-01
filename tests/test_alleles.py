import pandas as pd
import numpy as np
from pandas.testing import assert_frame_equal

from pgkb_hap_parser.parser import parse_alleles

def test_parse_alleles():
    """
    Input:
               0    1    2
    0  Reference    G    C
    1   c.520C>T    A  NaN
    2  c.3257G>A  NaN    T

    Output:
            name    alleles
    0  Reference      G,C
    1   c.520C>T      A,C
    2  c.3257G>A      G,T
    """
    obs = parse_alleles(
        pd.DataFrame(
            data=[
                ['Reference', 'G', 'C'],
                ['c.520C>T', 'A', np.NaN], 
                ['c.3257G>A', np.NaN, 'T']
            ]
        ),
        num_variants=2
     )
    exp = pd.DataFrame(
        data={
            'name': ['Reference', 'c.520C>T', 'c.3257G>A'],
            'alleles': ['G,C', 'A,C', 'G,T']
        }
    )
    assert_frame_equal(obs, exp)
    
    # with notes i.e

    # with function i.e G6PD
