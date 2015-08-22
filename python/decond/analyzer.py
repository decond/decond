import h5py
import numpy as np
import scipy.integrate as integrate
from ._version import __version__


class CorrFile(h5py.File):
    """
    Correlation data file output from decond.f90

    Can only open eithr a new file to write or an old file to read
    """
    def __init__(self, name, mode='r', **kwarg):
        if mode not in ('r', 'w-', 'x'):
            raise Error(type(self).__name__ +
                        " can only be opened in 'r', 'w-', 'x' mode")
        super().__init__(name, mode, **kwarg)
        self.filemode = mode
        self.buffer = CorrFile._Buffer()
        if 'spatialDec' in self:
            self.buffer.spatialDec = CorrFile._Buffer()
        else:
            self.buffer.spatialDec = None

        if 'energyDec' in self:
            self.buffer.energyDec = CorrFile._Buffer()
        else:
            self.buffer.energyDec = None

        if mode is 'r':
            self._check()
            self._read_corr_buffer()

    def _check(self):
        if 'version' in self.attrs:
            self.buffer.version = self.attrs['version']
            (fmajor, fminor, fpatch) = (self.buffer.version.decode().
                                        split(sep='.'))
            (major, minor, patch) = (__version__.split(sep='.'))
            if fmajor < major:
                raise Error(
                        "File " + self.filename +
                        " is of version " + fmajor + ".X.X, " +
                        "while this program requires at least " +
                        major + ".X.X")
        else:
            raise Error(
                    "File " + self.filename +
                    " has no version number. " +
                    "This program requires files of at least " +
                    "version " + major + ".X.X")

        filetype = self.attrs['type'].decode()
        if filetype != type(self).__name__:
            raise Error("Expecting {0} but {1} is encountered".format(
                type(self).__name__, filetype))

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

        if self.buffer.spatialDec is not None:
            self._read_corrdec_buffer('spatialDec')

        if self.buffer.energyDec is not None:
            self._read_corrdec_buffer('energyDec')

    def _read_corrdec_buffer(self, dectype):
        dec_group = self[dectype]
        buf = getattr(self.buffer, dectype)
        buf.decBins = dec_group['decBins'][...]
        buf.decBins_unit = dec_group['decBins'].attrs['unit']
        buf.decCorr = dec_group['decCorr'][...]
        buf.decCorr_unit = dec_group['decCorr'].attrs['unit']
        buf.decPairCount = dec_group['decPairCount'][...]

    def _cal_cesaro(self):
        """
        Calculate Cesaro data
        """
        def cesaro_integrate(y, x):
            cesaro = integrate.cumtrapz(y, x, initial=0)
            cesaro = integrate.cumtrapz(cesaro, x, initial=0)
            return cesaro

        self.buffer.nDCesaro = cesaro_integrate(self.buffer.nCorr,
                                                self.buffer.timeLags)
        self.buffer.nDCesaro_unit = self.buffer.nCorr_unit

        if self.buffer.spatialDec is not None:
            self.buffer.spatialDec.decDCesaro = cesaro_integrate(
                    self.buffer.spatialDec.decCorr, self.buffer.timeLags)
            self.buffer.spatialDec.decDCesaro_unit = (
                    self.buffer.spatialDec.decCorr_unit)

        if self.buffer.energyDec is not None:
            self.buffer.energyDec.decDCesaro = cesaro_integrate(
                    self.buffer.energyDec.decCorr, self.buffer.timeLags)
            self.buffer.energyDec.decDCesaro_unit = (
                    self.buffer.energyDec.decCorr_unit)

    class _Buffer():
        pass


class DecondFile(CorrFile):
    """
    Analyzed data
    """
    def __init__(self, name, mode='r', **kwarg):
        super().__init__(name, mode, **kwarg)
        if mode in ('r'):
            self._read_decond_buffer()
        else:
            self.buffer.numSample = 0  # initialize empty data

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
            self._read_deconddec_buffer('spatialDec')

        if 'energyDec' in self:
            self._read_deconddec_buffer('energyDec')

    def _read_deconddec_buffer(self, dectype):
        dec_group = self[dectype]
        buf = getattr(self.buffer, dectype)
        buf.decCorr_err = dec_group['decCorr_err'][...]
        buf.decPairCount_err = dec_group['decPairCount_err'][...]
        buf.decDCesaro = h5g_to_dict(dec_group['decDCesaro'])
        buf.decDCesaro_err = h5g_to_dict(dec_group['decDCesaro_err'])
        buf.decDCesaro_unit = self['decDCesaro'].attrs['unit']
        buf.decD = h5g_to_dict(dec_group['decD'])
        buf.decD_err = h5g_to_dict(dec_group['decD_err'])
        buf.decD_unit = self['decD'].attrs['unit']

    def _add_sample(self, samples):
        if not isinstance(samples, list):
            samples = [samples]

        if self.buffer.numSample == 0:
            with CorrFile(samples[0]) as cf:
                cf._cal_cesaro()
                self.buffer = cf.buffer
            self.buffer.numSample = 1

            def init_Err(buf, data_name):
                data_name_m2 = data_name + '_m2'
                data_name_err = data_name + '_err'
                setattr(buf, data_name_m2,
                        np.zeros_like(getattr(buf, data_name)))
                setattr(buf, data_name_err,
                        _m2_to_err(getattr(buf, data_name_m2), buf.numSample))

            init_Err(self.buffer, 'volume')
            init_Err(self.buffer, 'nCorr')
            init_Err(self.buffer, 'nDCesaro')

            def init_decErr(buf, num_sample):
                buf.decCorr_m2 = np.zeros_like(buf.decCorr)
                buf.decCorr_err = _m2_to_err(buf.decCorr_m2, num_sample)

                buf.decPairCount_m2 = np.zeros_like(buf.decPairCount)
                buf.decPairCount_err = _m2_to_err(
                        buf.decPairCount_m2, num_sample)

                buf.decDCesaro_m2 = np.zeros_like(buf.decDCesaro)
                buf.decDCesaro_err = _m2_to_err(buf.decDCesaro_m2, num_sample)

            if self.buffer.spatialDec is not None:
                init_decErr(self.buffer.spatialDec, self.buffer.numSample)

            if self.buffer.energyDec is not None:
                init_decErr(self.buffer.energyDec, self.buffer.numSample)

            begin = 1
        else:
            # not new file, buffer should have already been loaded
            self.buffer.volume_m2 = _err_to_m2(self.buffer.volume_err,
                                               self.buffer.numSample)
            self.buffer.nCorr_m2 = _err_to_m2(self.buffer.nCorr_err,
                                              self.buffer.numSample)
            self.buffer.nDCesaro_m2 = _err_to_m2(self.buffer.nDCesaro_err,
                                                 self.buffer.numSample)

            def init_dec_m2(buf, num_sample):
                buf.decCorr_m2 = _err_to_m2(buf.decCorr_err, num_sample)
                buf.decPairCount_m2 = _err_to_m2(buf.decPairCount_err,
                                                 num_sample)
                buf.decDCesaro_m2 = _err_to_m2(buf.decDCesaro_err, num_sample)

            if self.buffer.spatialDec is not None:
                init_dec_m2(self.buffer.spatialDec, self.buffer.numSample)

            if self.buffer.energyDec is not None:
                init_dec_m2(self.buffer.energyDec, self.buffer.numSample)

            begin = 0

        for sample in samples[begin:]:
            with CorrFile(sample) as f:
                self.buffer.numSample += 1
                f._cal_cesaro()
                self._add_data(self.buffer, self.buffer.numSample,
                               'volume', f.buffer.volume)
                self._add_data(self.buffer, self.buffer.numSample,
                               'nCorr', f.buffer.nCorr)
                self._add_data(self.buffer, self.buffer.numSample,
                               'nDCesaro', f.buffer.nDCesaro)

                def add_dec_data(buf, num_sample, new_buf):
                    self._add_data(buf, num_sample, 'decCorr', new_buf.decCorr)
                    self._add_data(buf, num_sample, 'decPairCount',
                                   new_buf.decPairCount)
                    self._add_data(buf, num_sample, 'decDCesaro',
                                   new_buf.decDCesaro)

                if self.buffer.spatialDec is not None:
                    add_dec_data(self.buffer.spatialDec, self.buffer.numSample,
                                 f.buffer.spatialDec)

                if self.buffer.energyDec is not None:
                    add_dec_data(self.buffer.energyDec, self.buffer.numSample,
                                 f.buffer.energyDec)

        self._fit_cesaro()

    @staticmethod
    def _add_data(buffer, num_sample, data_name, new_data):
        """
        Update the number, mean, and m2 of buffer.<data_name>,
        and add buffer.<data_name>_err

        http://www.wikiwand.com/en/Algorithms_for_calculating_variance#/On-line_algorithm
        """
        mean = getattr(buffer, data_name)
        m2 = getattr(buffer, data_name + '_m2')

        delta = new_data - mean
        mean += delta / num_sample
        m2 += delta * (new_data - mean)

        setattr(buffer, data_name + '_err',
                _m2_to_err(m2, num_sample))

    def _fit_cesaro(self):
        pass

    def close(self):
        if self.filemode in ('w-', 'x'):
            self._write_buffer()
        super().close()

    def _write_buffer(self):
        self.attrs['version'] = __version__
        self.attrs['type'] = np.string_(type(self).__name__)
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
#        self['nD'] = self.buffer.nD
#        self['nD_err'] = self.buffer.nD_err
#        self['nDCesaroTotal'] = self.buffer.nDCesaroTotal
#        self['nDCesaroTotal_err'] = self.buffer.nDCesaroTotal_err
#        self['nDTotal'] = self.buffer.nDTotal
#        self['nDTotal_err'] = self.buffer.nDTotal_err

        if self.buffer.spatialDec is not None:
            self._write_dec_buffer('spatialDec')

        if self.buffer.energyDec is not None:
            self._write_dec_buffer('energyDec')

    def _write_dec_buffer(self, dectype):
        dec_group = self.require_group(dectype)
        buf = getattr(self.buffer, dectype)
        dec_group['decBins'] = buf.decBins
        dec_group['decCorr'] = buf.decCorr
        dec_group['decPairCount'] = buf.decPairCount

        dec_group['decBins'].attrs['unit'] = buf.decBins_unit
        dec_group['decCorr'].attrs['unit'] = buf.decCorr_unit
        dec_group['decCorr_err'] = buf.decCorr_err
        dec_group['decPairCount_err'] = buf.decPairCount_err

        dec_group['decDCesaro'] = buf.decDCesaro
        dec_group['decDCesaro_err'] = buf.decDCesaro_err
        dec_group['decDCesaro'].attrs['unit'] = buf.decDCesaro_unit

#        dec_group['decD'] = buf.decD
#        dec_group['decD_err'] = buf.decD_err
#        dec_group['decD'].attrs['unit'] = buf.decD_unit


class Error(Exception):
    pass


def h5g_to_dict(group):
    D = {}
    for k, v in group.items():
        D[k] = v[...]
    return D


# def _get_filetype(file):
#     with h5py.File(file, 'r') as f:
#         return f.attrs['type'].decode()


def _err_to_m2(err, n):
    if n > 1:
        return np.square(err) * n * (n - 1)
    else:
        return np.zeros_like(err)


def _m2_to_err(m2, n):
    if n > 1:
        return np.sqrt(m2 / ((n - 1) * n))
    else:
        return np.full(m2.shape, np.nan)


def cal_decond(outname, samples):
    with DecondFile(outname, 'x') as outfile:
        outfile._add_sample(samples)
        return outfile.buffer
