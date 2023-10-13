# -*- coding: utf-8 -*-
""" Test for combineMultipleModels_1d
"""
from tempfile import NamedTemporaryFile
import os.path
import matplotlib as mpl
from matplotlib.testing.compare import compare_images
import baredSC.combineMultipleModels_1d

mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

EXAMPLE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..",  "..", "example")

BARED_1D_IMAGES_SUFFIX = ['', '_individuals',
                          '_with_posterior', '_posterior_individuals',
                          '_posterior_andco']

BARED_1D_TEXT_SUFFIX = ['_pdf.txt', '_posterior_per_cell.txt']

TOLERENCE = 13  # default matplotlib pixed difference tolerance


def test_combine_1d_test1():
    """Simple test from small npz
    """

    extension = 'png'
    expected = os.path.join(ROOT, 'combine_test1')
    with NamedTemporaryFile(suffix= f'.{extension}', prefix='baredSC_test_',
                            delete=False) as outfig:
        outfig_base = outfig.name[:-(len(extension) + 1)]
        args = f"--input {ROOT}/nih3t3_generated_2d_2.txt " \
            "--geneColName 0.5_0_0_0.5_x " \
            f"--outputs {ROOT}/small_1gauss {ROOT}/small_2gauss " \
            "--nx 10 --prettyBins 100 " \
            f"--figure {outfig.name}".split()
        args += ['--title', 'first gene combine 1 and 2 gauss']
        baredSC.combineMultipleModels_1d.main(args)
        for suffix in BARED_1D_IMAGES_SUFFIX:
            expected_file = f'{expected}{suffix}.{extension}'
            obtained_file = f'{outfig_base}{suffix}.{extension}'
            res = compare_images(expected_file,
                                obtained_file, TOLERENCE)
            assert res is None, res

            os.remove(obtained_file)
    # for suffix in BARED_1D_TEXT_SUFFIX:
    #     expected_file = f'{expected}{suffix}'
    #     obtained_file = f'{outfig_base}{suffix}'
    #     expected_mat = np.loadtxt(expected_file, skiprows=1)
    #     obtained_mat = np.loadtxt(obtained_file, skiprows=1)
    #     assert np.all(np.isclose(obtained_mat, expected_mat))

    #     os.remove(obtained_file)
