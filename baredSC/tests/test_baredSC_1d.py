# -*- coding: utf-8 -*-
""" Test for baredSC_1d
"""
from tempfile import NamedTemporaryFile
import os.path
import matplotlib as mpl
from matplotlib.testing.compare import compare_images
import baredSC.baredSC_1d

mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

EXAMPLE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..",  "..", "example")

BARED_1D_IMAGES_SUFFIX = ['', '_convergence', '_p', '_corner', '_individuals',
                          '_with_posterior', '_posterior_individuals',
                          '_posterior_andco']

BARED_1D_TEXT_SUFFIX = ['_neff.txt', '_p.txt', '_pdf.txt', '_posterior_per_cell.txt']

TOLERENCE = 18  # default matplotlib pixed difference tolerance


def test_baredSC_1d_1gauss_default():
    """First test matching example
    """
    extension = 'png'
    expected = os.path.join(EXAMPLE, 'first_example_1d_1gauss')
    outfile = NamedTemporaryFile(prefix='baredSC_test_',  # pylint: disable=R1732
                                 delete=False)
    outfig = NamedTemporaryFile(suffix= f'.{extension}', prefix='baredSC_test_',  # pylint: disable=R1732
                                delete=False)
    outfig_base = outfig.name[:-(len(extension) + 1)]
    outfile_evid = NamedTemporaryFile(suffix='.txt', prefix='baredSC_test_',  # pylint: disable=R1732
                                      delete=False)
    args = f"--input {ROOT}/nih3t3_generated_2d_2.txt " \
           "--geneColName 0.5_0_0_0.5_x " \
           f"--output {outfile.name}.npz " \
           "--nnorm 1 " \
           f"--figure {outfig.name} " \
           f"--logevidence {outfile_evid.name}".split()
    args += ['--title', 'first gene 1 gauss']
    baredSC.baredSC_1d.main(args)
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


def test_baredSC_1d_2gauss_log_pdf():
    """Second test with pdf cells subset prettyBins...
    """

    extension = 'pdf'
    expected = os.path.join(EXAMPLE, 'first_example_1d_2gauss_log')
    outfile = NamedTemporaryFile(prefix='baredSC_test_',  # pylint: disable=R1732
                                 delete=False)
    outfig = NamedTemporaryFile(suffix= f'.{extension}', prefix='baredSC_test_',  # pylint: disable=R1732
                                delete=False)
    outfig_base = outfig.name[:-(len(extension) + 1)]
    outfile_evid = NamedTemporaryFile(suffix='.txt', prefix='baredSC_test_',  # pylint: disable=R1732
                                      delete=False)
    args = f"--input {ROOT}/nih3t3_generated_2d_2.txt " \
           "--metadata1ColName 0_0.5_0.5_0_group " \
           "--metadata1Values 1.0 " \
           "--xmin -15 " \
           "--xmax -7 " \
           "--nx 25 " \
           "--xscale log " \
           "--minNeff 400 " \
           "--geneColName 0.5_0_0_0.5_x " \
           f"--output {outfile.name}.npz " \
           "--nnorm 2 " \
           f"--figure {outfig.name} " \
           f"--logevidence {outfile_evid.name}".split()
    args += ['--title', 'first gene 2 gauss log scale']
    baredSC.baredSC_1d.main(args)
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
