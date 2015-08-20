import h5py
import os.path
import numpy as np
import scipy as sp
from ._version import __version__


class CorrFile(h5py.File):
    """
    Correlation data file output from decond.f90
    """
    def __init__(self, name, mode='r', **kwarg):
        if isinstance(self, CorrFile) and (mode is not 'r'):
            raise IOError("CorrFile can only be opened in 'r' mode")
        file_exist_at_first = os.path.isfile(name)
        super().__init__(name, mode, **kwarg)
        self.buffer = _Buffer()
        if file_exist_at_first and mode in ('r', 'r+', 'a'):
            self._check_version()
            self._read_corr_buffer()

    def _check_version(self):
        if 'version' in self.attrs:
            self.buffer.version = self.attrs['version']
            (fmajor, fminor, fpatch) = (self.buffer.version.decode().
                                        split(sep='.'))
            (major, minor, patch) = (__version__.split(sep='.'))
            if fmajor < major:
                raise IOError(
                        "File " + self.filename +
                        " is of version " + fmajor + ".X.X, " +
                        "while this program requires at least " +
                        major + ".X.X")
        else:
            raise IOError(
                    "File " + self.filename +
                    " has no version number. " +
                    "This program requires files of at least " +
                    "version " + major + ".X.X")

    def _read_corr_buffer(self):
        self.buffer.charge = self['charge'][...]
        self.buffer.charge_unit = self['charge'].attrs['unit']
        self.buffer.numMol = self['numMol'][...]
        self.buffer.volume = self['volume'][...]
        self.buffer.volume_unit = self['volume'].attrs['unit']
        self.buffer.timeLags = self['timeLags'][...]
        self.buffer.timeLags_unit = self['timeLags'].attrs['unit']
        self.buffer.nCorr = self['nCorr'][...]
        self.buffer.nCorr_unit = self['nCorr'].attrs['unit']

        if 'spatialDec' in self:
            self.buffer.spatialDec = _Buffer()
            self._read_corrdec_buffer(self.buffer.spatialDec, 'spatialDec')

        if 'energyDec' in self:
            self.buffer.energyDec = _Buffer()
            self._read_corrdec_buffer(self.buffer.energyDec, 'spatialDec')

    def _read_corrdec_buffer(self, buf, dec_type):
        dec_group = self[dec_type]
        buf.decBins = dec_group['decBins'][...]
        buf.decBins_unit = dec_group['decBins'].attrs['unit']
        buf.decCorr = dec_group['decCorr'][...]
        buf.decCorr_unit = dec_group['decCorr'].attrs['unit']
        buf.decPairCount = dec_group['decPairCount'][...]

    def read_buffer(self):
        self._read_corr_buffer()


class DecondFile(CorrFile):
    """
    A container file of all the analyzed data
    """
    def __init__(self, name, mode='a', **kwarg):
        file_exist_at_first = os.path.isfile(name)
        super().__init__(name, mode, **kwarg)
        if file_exist_at_first and mode in ('r', 'r+', 'a'):
            self._read_decond_buffer()
        else:
            self._init_emptydata()

    def _init_emptydata(self):
        self.buffer.numSample = 0

    def _read_decond_buffer(self):
        self.buffer.numSample = self['numSample'][...]
        self.buffer.volume_err = self['volume_err'][...]
        self.buffer.nCorr_err = self['nCorr_err'][...]
        self.buffer.nDCesaro = h5g_to_dict(self['nDCesaro'])
        self.buffer.nDCesaro_err = h5g_to_dict(self['nDCesaro_err'])
        self.buffer.nDCesaro_unit = self['nDCesaro'].attrs['unit']
        self.buffer.nD = h5g_to_dict(self['nD'])
        self.buffer.nD_err = h5g_to_dict(self['nD_err'])
        self.buffer.nD_unit = self['nD'].attrs['unit']
        self.buffer.nDCesaroTotal = h5g_to_dict(self['nDCesaroTotal'])
        self.buffer.nDCesaroTotal_err = h5g_to_dict(
                self['nDCesaroTotal_err'])
        self.buffer.nDCesaroTotal_unit = self['nDCesaroTotal'].attrs['unit']
        self.buffer.nDTotal = h5g_to_dict(self['nDTotal'])
        self.buffer.nDTotal_err = h5g_to_dict(self['nDTotal_err'])
        self.buffer.nDTotal_unit = self['nDTotal'].attrs['unit']

        if 'spatialDec' in self:
            self._read_deconddec_buffer(self.buffer.spatialDec, 'spatialDec')

        if 'energyDec' in self:
            self._read_deconddec_buffer(self.buffer.energyDec, 'energyDec')

    def _read_deconddec_buffer(self, buf, dec_type):
        dec_group = self[dec_type]
        buf.decCorr_err = dec_group['decCorr_err'][...]
        buf.decPairCount_err = dec_group['decPairCount_err'][...]
        buf.decDCesaro = h5g_to_dict(dec_group['decDCesaro'])
        buf.decDCesaro_err = h5g_to_dict(dec_group['decDCesaro_err'])
        buf.decDCesaro_unit = self['decDCesaro'].attrs['unit']
        buf.decD = h5g_to_dict(dec_group['decD'])
        buf.decD_err = h5g_to_dict(dec_group['decD_err'])
        buf.decD_unit = self['decD'].attrs['unit']

    def read_buffer(self):
        self._read_corr_buffer()
        self._read_decond_buffer()

    def add_sample(self, samples):
        if not isinstance(samples, list):
            samples = [samples]

        if self.buffer.numSample == 0:
            # new file, initialize according to the first sample
            with CorrFile(samples[0]) as cf:
                self.buffer = _corr_to_cesaro_buffer(cf.buffer)
            self.buffer.numSample = 1
            begin = 1
        else:
            # read data from the existing decond file
            self.read_buffer()
            begin = 0

        for sample in samples[begin:]:
            with CorrFile(sample) as cf:
                pass

        self.write_buffer()

    def write_buffer(self):
        self.attrs['version'] = self.buffer.version
        self['numSample'] = self.buffer.numSample
        self['charge'] = self.buffer.charge
        self['numMol'] = self.buffer.numMol
        self['volume'] = self.buffer.volume
        self['volume'].attrs['unit'] = self.buffer.volume_unit
        self['volume_err'] = self.buffer.volume_err
        self['timeLags'] = self.buffer.timeLags
        self['nCorr'] = self.buffer.nCorr
        self['nCorr_err'] = self.buffer.nCorr_err
        self['nDCesaro'] = self.buffer.nDCesaro
        self['nDCesaro_err'] = self.buffer.nDCesaro_err
        self['nD'] = self.buffer.nD
        self['nD_err'] = self.buffer.nD_err
        self['nDCesaroTotal'] = self.buffer.nDCesaroTotal
        self['nDCesaroTotal_err'] = self.buffer.nDCesaroTotal_err
        self['nDTotal'] = self.buffer.nDTotal
        self['nDTotal_err'] = self.buffer.nDTotal_err

        if self.buffer.spatialDec is not None:
            self.buffer.spatialDec.write(self, 'spatialDec')

        if self.buffer.energyDec is not None:
            self.buffer.energyDec.write(self, 'energyDec')

    class DecData:
        def __init__(self, base, group):
            if group in base:
                dec_group = base['group']
            else:
                return None

            self.decBins = dec_group['decBins'][:]
            self.decCorr = dec_group['decCorr'][:, :, :]
            self.decPairCount = dec_group['decPairCount'][:, :]

            self.decBins_unit = DecondFile._readAttrStr(
                    dec_group['decBins'], 'unit')

            self.decCorr_unit = DecondFile._readAttrStr(
                    dec_group['decCorr'], 'unit')

            self.decCorr_err = DecondFile._readDset(dec_group, 'decCorr_err')

            self.decPairCount_err = DecondFile._readDset(
                    dec_group, 'decPairCount_err')

            self.decDCesaro = DecondFile._readGroupToDict(
                    dec_group, 'decDCesaro')
            self.decDCesaro_err = DecondFile._readGroupToDict(
                    dec_group, 'decDCesaro_rr')
            if self.decDCesaro is not None:
                self.decDCesaro_unit = DecondFile._readAttrStr(
                        dec_group['decDCesaro'], 'unit')

            self.decD = DecondFile._readGroupToDict(dec_group, 'decD')
            self.decD_err = DecondFile._readGroupToDict(dec_group, 'decD_err')

            if self.decD is not None:
                self.decD_unit = DecondFile._readAttrStr(
                        dec_group['decD'], 'unit')

        def write(self, base, group):
            dec_group = base.require_group(group)
            dec_group['decBins'] = self.decBins
            dec_group['decCorr'] = self.decCorr
            dec_group['decPairCount'] = self.decPairCount

            dec_group['decBins'].attrs['unit'] = self.decBins_unit
            dec_group['decCorr'].attrs['unit'] = self.decCorr_unit
            dec_group['decCorr_err'] = self.decCorr_err

            dec_group['decPairCount_err'] = self.decPairCount_err

            dec_group['decDCesaro'] = self.decDCesaro
            dec_group['decDCesaro_err'] = self.decDCesaro_err
            dec_group['decDCesaro'].attrs['unit'] = self.decDCesaro_unit

            dec_group['decD'] = self.decD
            dec_group['decD_err'] = self.decD_err
            dec_group['decD'].attrs['unit'] = self.decD_unit


class _Buffer():
    pass


class IOError(Exception):
    pass


def h5g_to_dict(group):
    D = {}
    for k, v in group.items():
        D[k] = v[...]
    return D


def _corr_to_cesaro_buffer(buf):
    """
    Return a buffer with Cesaro data added

    Input: a buffer with attributes
    buffer.timeLags (_unit)
    buffer.nCorr (_unit)
    [buffer.spatialDec.decCorr]
    [buffer.energyDec.decCorr]

    Output: a buffer with more attributes
    buffer.nDCesaro
    [buffer.spatialDec.decDCesaro]
    [buffer.energyDec.decDCesaro]
    """
    def cesaro_integrate(y, x):
        cesaro = sp.integrate.cumtrapz(y, x, initial=0)
        cesaro = sp.integrate.cumtrapz(cesaro, x, initial=0)
        return cesaro

    buf.nDCesaro = cesaro_integrate(buf.nCorr, buf.timeLags)
    buf.nDCesaro_unit = buf.nCorr_unit

    if hasattr(buf, 'spatialDec'):
        buf.spatialDec.decDCesaro = cesaro_integrate(
                buf.spatialDec.decCorr, buf.timeLags)
        buf.spatialDec.decDCesaro_unit = buf.spatialDec.decCorr_unit

    if hasattr(buf, 'energyDec'):
        buf.energyDec.decDCesaro = cesaro_integrate(
                buf.energyDec.decCorr, buf.timeLags)
        buf.energyDec.decDCesaro_unit = buf.energyDec.decCorr_unit


def err_to_m2(err, n):
    return np.square(err) * n * (n - 1)
