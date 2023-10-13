# -*- coding: utf-8 -*-
""" Test for combineMultipleModels_2d
"""
from tempfile import NamedTemporaryFile
import os.path
import matplotlib as mpl
from matplotlib.testing.compare import compare_images
import baredSC.combineMultipleModels_2d

mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

EXAMPLE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..",  "..", "example")

BARED_2D_IMAGES_SUFFIX = ['', '_individuals',
                          '_median']

BARED_2D_TEXT_SUFFIX = ['_corr.txt', '_pdf2d.txt', '_pdf2d_flat.txt']

TOLERENCE = 13  # default matplotlib pixed difference tolerance


def test_combine_2d_test1():
    """Simple test from small npz
    """

    extension = 'pdf'
    expected = os.path.join(ROOT, '2d_small_combined')
    with NamedTemporaryFile(suffix= f'.{extension}', prefix='baredSC_test_',
                            delete=False) as outfig:
        outfig_base = outfig.name[:-(len(extension) + 1)]
        args = f"--input {ROOT}/nih3t3_generated_2d_2.txt " \
           "--geneXColName 0.5_0_0_0.5_x " \
           "--geneYColName 0.5_0_0_0.5_y " \
           "--nx 10 --ny 12 " \
           "--prettyBinsx 50 --prettyBinsy 50 " \
            f"--outputs {ROOT}/2d_small_1gauss {ROOT}/2d_small_2gauss " \
            f"--figure {outfig.name}".split()
        baredSC.combineMultipleModels_2d.main(args)
        for suffix in BARED_2D_IMAGES_SUFFIX:
            expected_file = f'{expected}{suffix}.{extension}'
            obtained_file = f'{outfig_base}{suffix}.{extension}'
            res = compare_images(expected_file,
                                obtained_file, TOLERENCE)
            assert res is None, res

            os.remove(obtained_file)
    # for suffix in BARED_2D_TEXT_SUFFIX:
    #     expected_file = f'{expected}{suffix}'
    #     obtained_file = f'{outfig_base}{suffix}'
    #     expected_mat = np.loadtxt(expected_file, skiprows=1)
    #     obtained_mat = np.loadtxt(obtained_file, skiprows=1)
    #     assert np.all(np.isclose(obtained_mat, expected_mat))

    #     os.remove(obtained_file)
