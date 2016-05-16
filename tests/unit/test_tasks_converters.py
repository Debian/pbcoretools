
import subprocess
import tempfile
import unittest
import logging
import gzip
import os.path as op
import os

from pbcore.io import (FastaReader, FastqReader, openDataSet, HdfSubreadSet,
                       SubreadSet, ConsensusReadSet)
import pbcore.data
from pbcommand.testkit import PbTestApp
from pbcommand.utils import which

from pbcoretools import pbvalidate

from base import get_temp_file

log = logging.getLogger(__name__)

DATA = op.join(op.dirname(op.dirname(__file__)), "data")
BARCODED_SUBREAD_SET = op.join(DATA, "barcoded.subreadset.xml")


class Constants(object):
    BAX2BAM = "bax2bam"
    BAM2FASTA = "bam2fasta"
    BAM2BAM = "bam2bam"


SIV_DATA_DIR = "/pbi/dept/secondary/siv/testdata"


def _to_skip_msg(exe):
    return "Missing {e} or {d}".format(d=SIV_DATA_DIR, e=exe)

# XXX hacks to make sure tools are actually available
HAVE_BAX2BAM = which(Constants.BAX2BAM) is not None
HAVE_BAM2BAM = which(Constants.BAM2BAM) is not None
HAVE_BAM2FASTX = which(Constants.BAM2FASTA) is not None
HAVE_DATA_DIR = op.isdir(SIV_DATA_DIR)
HAVE_DATA_AND_BAX2BAM = HAVE_BAX2BAM and HAVE_DATA_DIR
HAVE_DATA_AND_BAM2BAM = HAVE_BAM2BAM and HAVE_DATA_DIR

SKIP_MSG_BAX2BAM = _to_skip_msg(Constants.BAX2BAM)
SKIP_MSG_BAM2FX = _to_skip_msg(Constants.BAM2FASTA)
SKIP_MSG_BAM2BAM = _to_skip_msg(Constants.BAM2BAM)

skip_unless_bax2bam = unittest.skipUnless(HAVE_DATA_AND_BAX2BAM, SKIP_MSG_BAX2BAM)
skip_unless_bam2fastx = unittest.skipUnless(HAVE_BAM2FASTX, SKIP_MSG_BAM2FX)
skip_unless_bam2bam = unittest.skipUnless(HAVE_DATA_AND_BAM2BAM, SKIP_MSG_BAM2BAM)


def _get_bax2bam_inputs():
    """Little hackery to get the setup class Inputs and to avoid calls to
    setupclass if skiptest is used

    Nat: we want to test that this behaves properly when multiple movies are
    supplied as input, so we make an HdfSubreadSet on the fly from various
    bax files in testdata
    """
    if HAVE_DATA_AND_BAX2BAM:
        hdf_subread_xml = tempfile.NamedTemporaryFile(suffix=".hdfsubreadset.xml").name

        bax_files = (SIV_DATA_DIR + "/SA3-RS/lambda/2372215/0007_tiny/Analysis_Results/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.bax.h5",
                     pbcore.data.getBaxH5_v23()[0])

        ds = HdfSubreadSet(*bax_files)
        ds.name = "lambda_rsii"
        assert len(set([f.movieName for f in ds.resourceReaders()])) == 2
        ds.write(hdf_subread_xml)
        return [hdf_subread_xml]
    else:
        # Assume the test data isn't found and the test won't be run
        return ["/path/to/this-test-should-be-skipped.txt"]


@skip_unless_bax2bam
class TestBax2Bam(PbTestApp):
    TASK_ID = "pbcoretools.tasks.h5_subreads_to_subread"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '

    # See comments above
    INPUT_FILES = _get_bax2bam_inputs()
    MAX_NPROC = 24

    RESOLVED_NPROC = 1
    RESOLVED_TASK_OPTIONS = {}
    IS_DISTRIBUTED = True
    RESOLVED_IS_DISTRIBUTED = True

    def run_after(self, rtc, output_dir):
        with SubreadSet(rtc.task.output_files[0]) as ds_out:
            self.assertEqual(len(ds_out.toExternalFiles()), 2)
            self.assertEqual(ds_out.name, "lambda_rsii")


@skip_unless_bam2bam
class TestBam2Bam(PbTestApp):
    TASK_ID = "pbcoretools.tasks.bam2bam_barcode"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [
        "/pbi/dept/secondary/siv/testdata/SA3-Sequel/phi29/315/3150101/r54008_20160219_002905/1_A01_micro/m54008_160219_003234_micro.subreadset.xml",
        "/pbi/dept/secondary/siv/barcodes/pacbio_barcodes_384/pacbio_barcodes_384.barcodeset.xml"
    ]
    MAX_NPROC = 8
    RESOLVED_NPROC = 8
    IS_DISTRIBUTED = True
    RESOLVED_IS_DISTRIBUTED = True

    def run_after(self, rtc, output_dir):
        err, metrics = pbvalidate.validate_dataset(
            file_name=rtc.task.output_files[0],
            dataset_type="SubreadSet",
            quick=False,
            validate_index=True,
            strict=True)
        self.assertEqual(len(err), 0)
        with SubreadSet(rtc.task.output_files[0]) as ds:
            self.assertEqual(len(ds.externalResources), 1)
            # make sure metadata are propagated
            md = ds.metadata
            self.assertEqual(
                md.collections.submetadata[0].attrib['InstrumentName'],
                "Inst54008")
            self.assertTrue(ds.externalResources[0].scraps is not None)
            self.assertEqual(ds.externalResources[0].barcodes,
                             self.INPUT_FILES[1])
            rr = ds.resourceReaders()[0]
            self.assertTrue(rr.pbi.hasBarcodeInfo)
            #self.assertEqual(len(rr.pbi.bcReverse), 13194)


class _BaseTestBam2Fasta(PbTestApp):
    TASK_ID = "pbcoretools.tasks.bam2fasta"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [get_temp_file(suffix=".subreadset.xml")]
    MAX_NPROC = 24
    RESOLVED_NPROC = 1
    IS_DISTRIBUTED = True
    RESOLVED_IS_DISTRIBUTED = True
    READER_CLASS = FastaReader

    def _get_output_file(self, rtc):
        return rtc.task.output_files[0]

    def _get_counts(self, rtc):
        with openDataSet(self.INPUT_FILES[0]) as ds:
            n_expected = len([rec for rec in ds])
        with self.READER_CLASS(self._get_output_file(rtc)) as f:
            n_actual = len([rec for rec in f])
        return n_expected, n_actual

    def run_after(self, rtc, output_dir):
        n_expected, n_actual = self._get_counts(rtc)
        self.assertEqual(n_actual, n_expected)


@skip_unless_bam2fastx
class TestBam2Fasta(_BaseTestBam2Fasta):

    @classmethod
    def setUpClass(cls):
        ds = SubreadSet(pbcore.data.getUnalignedBam(), strict=True)
        ds.write(cls.INPUT_FILES[0])
        super(TestBam2Fasta, cls).setUpClass()


@skip_unless_bam2fastx
class TestBam2FastaIgnoreBarcodes(_BaseTestBam2Fasta):
    """
    Make sure the base bam2fasta task always outputs a single FASTA file
    even when barcoding is present.
    """
    INPUT_FILES = [BARCODED_SUBREAD_SET]


@skip_unless_bam2fastx
class TestBam2FastaNoFilter(TestBam2Fasta):
    TASK_ID = "pbcoretools.tasks.bam2fasta_nofilter"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)


@skip_unless_bam2fastx
class TestBam2Fastq(TestBam2Fasta):
    TASK_ID = "pbcoretools.tasks.bam2fastq"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    READER_CLASS = FastqReader


@skip_unless_bam2fastx
class TestBam2FastqFiltered(TestBam2Fastq):
    TASK_OPTIONS = {"pbcoretools.task_options.min_subread_length": 1000}
    RESOLVED_TASK_OPTIONS = {"pbcoretools.task_options.min_subread_length": 1000}

    def run_after(self, rtc, output_dir):
        n_expected, n_actual = self._get_counts(rtc)
        self.assertTrue(0 < n_actual < n_expected,
            "FAILED: 0 < {a} < {e}".format(a=n_actual, e=n_expected))


@skip_unless_bam2fastx
class TestBam2FastaArchive(TestBam2Fasta):
    TASK_ID = "pbcoretools.tasks.bam2fasta_archive"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    TASK_OPTIONS = {"pbcoretools.task_options.min_subread_length": 1000}
    RESOLVED_TASK_OPTIONS = {"pbcoretools.task_options.min_subread_length": 1000}

    def _get_output_file(self, rtc):
        return gzip.open(rtc.task.output_files[0])

    def run_after(self, rtc, output_dir):
        n_expected, n_actual = self._get_counts(rtc)
        self.assertTrue(0 < n_actual < n_expected,
            "FAILED: 0 < {a} < {e}".format(a=n_actual, e=n_expected))


@skip_unless_bam2fastx
class TestBam2FastqArchive(TestBam2Fastq):
    TASK_ID = "pbcoretools.tasks.bam2fastq_archive"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)

    def _get_output_file(self, rtc):
        return gzip.open(rtc.task.output_files[0])


@skip_unless_bam2fastx
class TestBam2FastaCCS(TestBam2FastqArchive):
    TASK_ID = "pbcoretools.tasks.bam2fasta_ccs"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    INPUT_FILES = [get_temp_file(".consensusreadset.xml")]
    READER_CLASS = FastaReader

    @classmethod
    def setUpClass(cls):
        ds = ConsensusReadSet(pbcore.data.getCCSBAM(), strict=True)
        ds.write(cls.INPUT_FILES[0])


@skip_unless_bam2fastx
class TestBam2FastqCCS(TestBam2FastaCCS):
    TASK_ID = "pbcoretools.tasks.bam2fastq_ccs"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    READER_CLASS = FastqReader


@skip_unless_bam2fastx
class TestBam2FastaBarcoded(PbTestApp):
    TASK_ID = "pbcoretools.tasks.bam2fasta_archive"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [BARCODED_SUBREAD_SET]
    MAX_NPROC = 24
    RESOLVED_NPROC = 1
    IS_DISTRIBUTED = True
    RESOLVED_IS_DISTRIBUTED = True
    READER_CLASS = FastaReader
    EXT = "fasta"

    def run_after(self, rtc, output_dir):
        tmp_dir = tempfile.mkdtemp()
        _cwd = os.getcwd()
        try:
            os.chdir(tmp_dir)
            args = ["tar", "xzf", rtc.task.output_files[0]]
            self.assertEqual(subprocess.call(args), 0)
            file_names = sorted(os.listdir(tmp_dir))
            self.assertEqual(file_names,
                             ["reads.{e}.0_0.{e}".format(e=self.EXT),
                              "reads.{e}.2_2.{e}".format(e=self.EXT)])
            fastx_ids = ["m54008_160219_003234/74056024/3985_5421", # bc 0
                         "m54008_160219_003234/28901719/5058_5262" ] # bc 2
            for file_name, fastx_id in zip(file_names, fastx_ids):
                with self.READER_CLASS(file_name) as f:
                    records = [rec.id for rec in f]
                    self.assertEqual(len(records), 1)
                    self.assertEqual(records[0], fastx_id)
        finally:
            os.chdir(_cwd)


@skip_unless_bam2fastx
class TestBam2FastqBarcoded(TestBam2FastaBarcoded):
    TASK_ID = "pbcoretools.tasks.bam2fastq_archive"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    READER_CLASS = FastqReader
    EXT = "fastq"



class TestFasta2Fofn(PbTestApp):
    TASK_ID = "pbcoretools.tasks.fasta2fofn"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [pbcore.data.getLambdaFasta()]
    IS_DISTRIBUTED = False
    RESOLVED_IS_DISTRIBUTED = False


class TestFasta2ReferenceSet(PbTestApp):
    TASK_ID = "pbcoretools.tasks.fasta2referenceset"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [pbcore.data.getLambdaFasta()]
