"""
Functions for data conversion
==============================

"""

def py3_FixNPLoad(np):

    """
    Function to fix numpy.load in Python 3 environment
            -- requires allow_pickle = True
    """

    import numpy as np
    from functools import partial

    # save np.load
    np_load_old = partial(np.load)

    # modify the default parameters of np.load
    np.load = lambda *a,**k: np_load_old(*a, allow_pickle=True, encoding = 'latin1', **k)

    # call load_data with allow_pickle implicitly set to true
    # (train_data, train_labels), (test_data, test_labels) = imdb.load_data(num_words=10000)

    # restore np.load for future normal usage
    # np.load = np_load_old

    return np
