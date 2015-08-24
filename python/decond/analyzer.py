import h5py
import numpy as np
import scipy.integrate as integrate
from ._version import __version__


class DecType:
    spatial = 'spatialDec'
    energy = 'energyDec'


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
        if DecType.spatial in self:
            self.buffer.spatialDec = CorrFile._Buffer()
        else:
            self.buffer.spatialDec = None

        if DecType.energy in self:
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
        self.buffer.timeLags_width = (self.buffer.timeLags[1] -
                                      self.buffer.timeLags[0])
        self.buffer.nCorr = self['nCorr'][...]
        self.buffer.nCorr_unit = self['nCorr'].attrs['unit']

        def do_dec(dectype):
            dec_group = self[dectype]
            buf = getattr(self.buffer, dectype)
            buf.decBins = dec_group['decBins'][...]
            buf.decBins_unit = dec_group['decBins'].attrs['unit']
            buf.decBins_width = buf.decBins[1] - buf.decBins[0]
            buf.decCorr = dec_group['decCorr'][...]
            buf.decCorr_unit = dec_group['decCorr'].attrs['unit']
            buf.decPairCount = dec_group['decPairCount'][...]

        if self.buffer.spatialDec is not None:
            do_dec(DecType.spatial)

        if self.buffer.energyDec is not None:
            do_dec(DecType.energy)

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

        def do_dec(dectype):
            buf = getattr(self.buffer, dectype)
            buf.decDCesaro = cesaro_integrate(
                    buf.decCorr, self.buffer.timeLags)
            buf.decDCesaro_unit = (buf.decCorr_unit)

        if self.buffer.spatialDec is not None:
            do_dec(DecType.spatial)

        if self.buffer.energyDec is not None:
            do_dec(DecType.energy)

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

        def do_dec(dectype):
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

        if DecType.spatial in self:
            do_dec(DecType.spatial)

        if DecType.energy in self:
            do_dec(DecType.energy)

    def _add_sample(self, samples):
        if not isinstance(samples, list):
            samples = [samples]

        if self.buffer.numSample == 0:
            with CorrFile(samples[0]) as cf:
                cf._cal_cesaro()
                self.buffer = cf.buffer
            self.buffer.numSample = 1

            def init_Err(data_name):
                data_name_m2 = data_name + '_m2'
                data_name_err = data_name + '_err'
                setattr(self.buffer, data_name_m2,
                        np.zeros_like(getattr(self.buffer, data_name)))
                setattr(self.buffer, data_name_err,
                        _m2_to_err(getattr(self.buffer, data_name_m2),
                                   self.buffer.numSample))

            init_Err('volume')
            init_Err('nCorr')
            init_Err('nDCesaro')

            def init_decErr(dectype):
                buf = getattr(self.buffer, dectype)
                num_sample = self.buffer.numSample

                buf.decCorr_m2 = np.zeros_like(buf.decCorr)
                buf.decCorr_err = _m2_to_err(buf.decCorr_m2, num_sample)

                buf.decPairCount_m2 = np.zeros_like(buf.decPairCount)
                buf.decPairCount_err = _m2_to_err(
                        buf.decPairCount_m2, num_sample)

                buf.decDCesaro_m2 = np.zeros_like(buf.decDCesaro)
                buf.decDCesaro_err = _m2_to_err(buf.decDCesaro_m2, num_sample)

            if self.buffer.spatialDec is not None:
                init_decErr(DecType.spatial)

            if self.buffer.energyDec is not None:
                init_decErr(DecType.energy)

            begin = 1
        else:
            # not new file, buffer should have already been loaded
            self.buffer.volume_m2 = _err_to_m2(self.buffer.volume_err,
                                               self.buffer.numSample)
            self.buffer.nCorr_m2 = _err_to_m2(self.buffer.nCorr_err,
                                              self.buffer.numSample)
            self.buffer.nDCesaro_m2 = _err_to_m2(self.buffer.nDCesaro_err,
                                                 self.buffer.numSample)

            def init_dec_m2(dectype):
                buf = getattr(self.buffer, dectype)
                num_sample = self.buffer.numSample
                buf.decCorr_m2 = _err_to_m2(buf.decCorr_err, num_sample)
                buf.decPairCount_m2 = _err_to_m2(buf.decPairCount_err,
                                                 num_sample)
                buf.decDCesaro_m2 = _err_to_m2(buf.decDCesaro_err, num_sample)

            if self.buffer.spatialDec is not None:
                init_dec_m2(DecType.spatial)

            if self.buffer.energyDec is not None:
                init_dec_m2(DecType.energy)

            begin = 0

        # add more samples one by one
        for sample in samples[begin:]:
            with CorrFile(sample) as f:
                self.buffer.numSample += 1
                f._cal_cesaro()
                self._add_data('volume', f.buffer.volume)
                self._add_data('nCorr', f.buffer.nCorr)
                self._add_data('nDCesaro', f.buffer.nDCesaro)

                def add_dec_data(dectype, new_buf):
                    buf = getattr(new_buf, dectype)
                    self._add_data('decCorr', buf.decCorr, dectype=dectype)
                    self._add_data('decPairCount', buf.decPairCount,
                                   dectype=dectype)
                    self._add_data('decDCesaro', buf.decDCesaro,
                                   dectype=dectype)

                if self.buffer.spatialDec is not None:
                    add_dec_data(DecType.spatial, f.buffer)

                if self.buffer.energyDec is not None:
                    add_dec_data(DecType.energy, f.buffer)

        self._fit_cesaro()

    def _add_data(self, data_name, new_data, dectype=None):
        """
        Update the number, mean, and m2 of buffer.<data_name>,
        and add buffer.<data_name>_err

        http://www.wikiwand.com/en/Algorithms_for_calculating_variance#/On-line_algorithm
        """
        if dectype is None:
            buf = self.buffer
        else:
            buf = getattr(self.buffer, dectype)

        num_sample = self.buffer.numSample
        mean = getattr(buf, data_name)
        m2 = getattr(buf, data_name + '_m2')

        delta = new_data - mean
        mean += delta / num_sample
        m2 += delta * (new_data - mean)

        if np.isscalar(mean):
            setattr(buf, data_name, mean)
            setattr(buf, data_name + '_m2', m2)

        setattr(buf, data_name + '_err',
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

        def do_dec(dectype):
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

#            dec_group['decD'] = buf.decD
#            dec_group['decD_err'] = buf.decD_err
#            dec_group['decD'].attrs['unit'] = buf.decD_unit

        if self.buffer.spatialDec is not None:
            do_dec(DecType.spatial)

        if self.buffer.energyDec is not None:
            do_dec(DecType.energy)


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


def _get_inner_index(a, b):
    awidth = a[1] - a[0]
    bwidth = b[1] - b[0]
    if not np.isclose(awidth, bwidth):
        raise Error("Bin widths should be the same (within tolerance)")

    if not np.isclose(np.mod(a[0] - b[0], awidth), 0.0):
        raise Error("Both bins should have a common grid")

    (a_begin, a_end, b_begin, b_end) = (0, len(a) - 1, 0, len(b) - 1)

    if b[0] > a[0]:
        a_begin = np.round((b[0] - a[0]) / awidth)
    elif b[0] < a[0]:
        b_begin = np.round((a[0] - b[0]) / bwidth)

    if b[-1] < a[-1]:
        a_end = np.round((b[-1] - a[0]) / awidth)
    elif b[-1] > a[-1]:
        b_end = np.round((a[-1] - b[0]) / bwidth)

    return (a_begin, a_end), (b_begin, b_end)


def cal_decond(outname, samples):
    with DecondFile(outname, 'x') as outfile:
        outfile._add_sample(samples)
        return outfile.buffer
