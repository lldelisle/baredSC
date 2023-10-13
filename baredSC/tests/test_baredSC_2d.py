# -*- coding: utf-8 -*-
""" Test for baredSC_2d
"""
from tempfile import NamedTemporaryFile
import os.path
import matplotlib as mpl
from matplotlib.testing.compare import compare_images
import baredSC.baredSC_2d

mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

EXAMPLE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..",  "..", "example")

BARED_2D_IMAGES_SUFFIX = ['', '_convergence', '_p', '_corner', '_individuals',
                          '_median']

BARED_2D_TEXT_SUFFIX = ['_neff.txt', '_p.txt', '_corr.txt', '_pdf2d.txt', '_pdf2d_flat.txt']

TOLERENCE = 13  # default matplotlib pixed difference tolerance


def test_baredSC_2d_1gauss_small():
    """First test with small
    """
    extension = 'png'
    expected = os.path.join(ROOT, '2d_small_2gauss')
    outfile = NamedTemporaryFile(prefix='baredSC_test_',  # pylint: disable=R1732
                                 delete=False)
    outfig = NamedTemporaryFile(suffix= f'.{extension}', prefix='baredSC_test_',  # pylint: disable=R1732
                                delete=False)
    outfig_base = outfig.name[:-(len(extension) + 1)]
    outfile_evid = NamedTemporaryFile(suffix='.txt', prefix='baredSC_test_',  # pylint: disable=R1732
                                      delete=False)
    args = f"--input {ROOT}/nih3t3_generated_2d_2.txt " \
           "--geneXColName 0.5_0_0_0.5_x " \
           "--geneYColName 0.5_0_0_0.5_y " \
           "--nx 10 --ny 12 --nsampMCMC 20000 " \
           "--prettyBinsx 50 --prettyBinsy 50 " \
           f"--output {outfile.name}.npz " \
           "--nnorm 2 " \
           f"--figure {outfig.name} " \
           f"--logevidence {outfile_evid.name}".split()
    baredSC.baredSC_2d.main(args)
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
