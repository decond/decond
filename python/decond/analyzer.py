import h5py
import numpy as np
import scipy.integrate as integrate
from enum import Enum
from ._version import __version__


class DecType(Enum):
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

        if mode is 'r':
            for type_ in DecType:
                if type_.value in self:
                    setattr(self.buffer, type_.value, CorrFile._Buffer())
                else:
                    setattr(self.buffer, type_.value, None)

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
            dec_group = self[dectype.value]
            buf = getattr(self.buffer, dectype.value)
            buf.decBins = dec_group['decBins'][...]
            buf.decBins_unit = dec_group['decBins'].attrs['unit']
            buf.decBins_width = buf.decBins[1] - buf.decBins[0]
            buf.decCorr = dec_group['decCorr'][...]
            buf.decCorr_unit = dec_group['decCorr'].attrs['unit']
            buf.decPairCount = dec_group['decPairCount'][...]

        for type_ in DecType:
            if type_.value in self:
                do_dec(type_)

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

        # Unit: nCorr (L^2 T^-2), nDCesaro (L^2)
        self.buffer.nDCesaro_unit = np.string_(
                self.buffer.nCorr_unit.decode().split()[0])

        def do_dec(buf):
            buf.decDCesaro = cesaro_integrate(
                    buf.decCorr, self.buffer.timeLags)
            buf.decDCesaro_unit = np.string_(
                    buf.decCorr_unit.decode().split()[0])

        for type_ in DecType:
            buf = getattr(self.buffer, type_.value)
            if buf is not None:
                do_dec(buf)

    def _write_buffer(self):
        self.attrs['version'] = np.string_(__version__)
        self.attrs['type'] = np.string_(type(self).__name__)
        self['charge'] = self.buffer.charge
        self['charge'].attrs['unit'] = self.buffer.charge_unit
        self['numMol'] = self.buffer.numMol
        self['volume'] = self.buffer.volume
        self['volume'].attrs['unit'] = self.buffer.volume_unit
        self['timeLags'] = self.buffer.timeLags
        self['timeLags'].attrs['unit'] = self.buffer.timeLags_unit
        self['nCorr'] = self.buffer.nCorr
        self['nCorr'].attrs['unit'] = self.buffer.nCorr_unit

        def do_dec(dectype):
            dec_group = self.require_group(dectype.value)
            buf = getattr(self.buffer, dectype.value)
            dec_group['decBins'] = buf.decBins
            dec_group['decCorr'] = buf.decCorr
            dec_group['decPairCount'] = buf.decPairCount

            dec_group['decBins'].attrs['unit'] = buf.decBins_unit
            dec_group['decCorr'].attrs['unit'] = buf.decCorr_unit

        for type_ in DecType:
            if getattr(self.buffer, type_.value) is not None:
                do_dec(type_)

    def _shrink_corr_buffer(self, sel):
        self.buffer.timeLags = self.buffer.timeLags[sel]
        self.buffer.nCorr = self.buffer.nCorr[..., sel]

    def _shrink_dec_buffer(self, dectype, sel, sel_dec):
        buf = getattr(self.buffer, dectype.value)
        buf.decBins = buf.decBins[sel_dec]
        buf.decCorr = buf.decCorr[:, sel_dec, sel]
        buf.decPairCount = buf.decPairCount[:, sel_dec]

    def _intersect_buffer(self, new_file):
        s_sel, n_sel = _get_inner_index(
                self.buffer.timeLags, new_file.buffer.timeLags)

        self._shrink_corr_buffer(s_sel)
        new_file._shrink_corr_buffer(n_sel)
#        assert(np.allclose(self.buffer.timeLags, new_file.buffer.timeLags))

        def do_dec(dectype):
            sb_dec = getattr(self.buffer, dectype.value)
            nb_dec = getattr(new_file.buffer, dectype.value)
            s_sel_dec, n_sel_dec = _get_inner_index(
                sb_dec.decBins, nb_dec.decBins)

            self._shrink_dec_buffer(dectype, s_sel, s_sel_dec)
            new_file._shrink_dec_buffer(dectype, n_sel, n_sel_dec)
#            assert(np.allclose(sb_dec.decBins, nb_dec.decBins))

        for type_ in DecType:
            if getattr(self.buffer, type_.value) is not None:
                do_dec(type_)

    def close(self):
        if self.filemode in ('w-', 'x'):
            self._write_buffer()
        super().close()

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
        self.buffer.nDCesaro = self['nDCesaro'][...]
        self.buffer.nDCesaro_err = self['nDCesaro_err'][...]
        self.buffer.nDCesaro_unit = self['nDCesaro'].attrs['unit']
        self.buffer.fit = self['fit'][...]
        self.buffer.fit_unit = self['fit'].attrs['unit']
#        self.buffer.nD = h5g_to_dict(self['nD'])
#        self.buffer.nD_err = h5g_to_dict(self['nD_err'])
#        self.buffer.nD_unit = self['nD'].attrs['unit']
#        self.buffer.nDCesaroTotal = self['nDCesaroTotal'][...]
#        self.buffer.nDCesaroTotal_err = self['nDCesaroTotal_err'][...]
#        self.buffer.nDCesaroTotal_unit = self['nDCesaroTotal'].attrs['unit']
#        self.buffer.nDTotal = h5g_to_dict(self['nDTotal'])
#        self.buffer.nDTotal_err = h5g_to_dict(self['nDTotal_err'])
#        self.buffer.nDTotal_unit = self['nDTotal'].attrs['unit']

        def do_dec(dectype):
            dec_group = self[dectype.value]
            buf = getattr(self.buffer, dectype.value)
            buf.decCorr_err = dec_group['decCorr_err'][...]
            buf.decPairCount_err = dec_group['decPairCount_err'][...]
            buf.decDCesaro = dec_group['decDCesaro'][...]
            buf.decDCesaro_err = dec_group['decDCesaro_err'][...]
            buf.decDCesaro_unit = dec_group['decDCesaro'].attrs['unit']
#            buf.decD = h5g_to_dict(dec_group['decD'])
#            buf.decD_err = h5g_to_dict(dec_group['decD_err'])
#            buf.decD_unit = self['decD'].attrs['unit']

        for type_ in DecType:
            if type_.value in self:
                do_dec(type_)

    def _add_sample(self, samples, fit):
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

            def init_decErr(buf):
                num_sample = self.buffer.numSample

                buf.decCorr_m2 = np.zeros_like(buf.decCorr)
                buf.decCorr_err = _m2_to_err(buf.decCorr_m2, num_sample)  # nan

                buf.decPairCount_m2 = np.zeros_like(buf.decPairCount)
                buf.decPairCount_err = _m2_to_err(
                        buf.decPairCount_m2, num_sample)  # nan

                buf.decDCesaro_m2 = np.zeros_like(buf.decDCesaro)
                buf.decDCesaro_err = _m2_to_err(buf.decDCesaro_m2, num_sample)

            for type_ in DecType:
                buf = getattr(self.buffer, type_.value)
                if buf is not None:
                    init_decErr(buf)

            begin = 1
        else:
            # not new file, buffer should have already been loaded
            self.buffer.volume_m2 = _err_to_m2(self.buffer.volume_err,
                                               self.buffer.numSample)
            self.buffer.nCorr_m2 = _err_to_m2(self.buffer.nCorr_err,
                                              self.buffer.numSample)
            self.buffer.nDCesaro_m2 = _err_to_m2(self.buffer.nDCesaro_err,
                                                 self.buffer.numSample)

            def init_dec_m2(buf):
                num_sample = self.buffer.numSample
                buf.decCorr_m2 = _err_to_m2(
                        buf.decCorr_err, num_sample, buf.decPairCount)
                buf.decDCesaro_m2 = _err_to_m2(
                        buf.decDCesaro_err, num_sample, buf.decPairCount)
                buf.decPairCount_m2 = _err_to_m2(
                        buf.decPairCount_err, num_sample)

            for type_ in DecType:
                buf = getattr(self.buffer, type_.value)
                if buf is not None:
                    init_dec_m2(buf)

            begin = 0

        # add more samples one by one
        def add_data(data_name, new_data, dectype=None):
            """
            Update the number, mean, and m2 of buffer.<data_name>,
            and add buffer.<data_name>_err

            http://www.wikiwand.com/en/Algorithms_for_calculating_variance#/On-line_algorithm
            """
            if dectype is None:
                buf = self.buffer
            else:
                buf = getattr(self.buffer, dectype.value)

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

        def add_weighted_data(data_name, weight_name, new_data, new_weight,
                              dectype=None):
            """
            http://www.wikiwand.com/en/Algorithms_for_calculating_variance#/Weighted_incremental_algorithm
            """
            if dectype is None:
                buf = self.buffer
            else:
                buf = getattr(self.buffer, dectype.value)

            old_weight = getattr(buf, weight_name)
            sum_weight = (self.buffer.numSample - 1) * old_weight
            mean = getattr(buf, data_name)
            m2 = getattr(buf, data_name + '_m2')

            temp = new_weight + sum_weight
            delta = new_data - mean
            r = delta * new_weight[..., np.newaxis] / temp[..., np.newaxis]
            mean += r
            m2 += sum_weight[..., np.newaxis] * delta * r
            sum_weight = temp

            if np.isscalar(mean):
                setattr(buf, data_name, mean)
                setattr(buf, data_name + '_m2', m2)

            # it is important to explicitly extend the sum_weight dimension
            # by [..., np.newaxis], otherwis, numpy will broadcast it
            # wrongly when the dimension of sum_weight is [N, N] (square)
            setattr(buf, data_name + '_err', _m2_to_err(
                m2, self.buffer.numSample,
                sum_weight[..., np.newaxis] / self.buffer.numSample))

        def add_dec_data(dectype, new_buf):
            buf = getattr(new_buf, dectype.value)
            add_weighted_data('decCorr', 'decPairCount',
                              buf.decCorr, buf.decPairCount, dectype)
            add_weighted_data('decDCesaro', 'decPairCount',
                              buf.decDCesaro, buf.decPairCount, dectype)

            # Note that decPairCount must be updated last
            add_data('decPairCount', buf.decPairCount, dectype)

        for sample in samples[begin:]:
            with CorrFile(sample) as f:
                self.buffer.numSample += 1
                self._intersect_buffer(f)
                f._cal_cesaro()

                add_data('volume', f.buffer.volume)
                add_data('nCorr', f.buffer.nCorr)
                add_data('nDCesaro', f.buffer.nDCesaro)

                for type_ in DecType:
                    if getattr(self.buffer, type_.value) is not None:
                        add_dec_data(type_, f.buffer)

        self._fit_cesaro(fit)

    def _shrink_corr_buffer(self, sel):
        super()._shrink_corr_buffer(sel)

        self.buffer.nCorr_m2 = self.buffer.nCorr_m2[..., sel]
        self.buffer.nCorr_err = self.buffer.nCorr_err[..., sel]
        self.buffer.nDCesaro = self.buffer.nDCesaro[..., sel]
        self.buffer.nDCesaro_m2 = self.buffer.nDCesaro_m2[..., sel]
        self.buffer.nDCesaro_err = self.buffer.nDCesaro_err[..., sel]

    def _shrink_dec_buffer(self, dectype, sel, sel_dec):
        super()._shrink_dec_buffer(dectype, sel, sel_dec)

        buf = getattr(self.buffer, dectype.value)
        buf.decCorr_m2 = buf.decCorr_m2[:, sel_dec, sel]
        buf.decCorr_err = buf.decCorr_err[:, sel_dec, sel]
        buf.decPairCount_m2 = buf.decPairCount_m2[:, sel_dec]
        buf.decPairCount_err = buf.decPairCount_err[:, sel_dec]
        buf.decDCesaro = buf.decDCesaro[:, sel_dec, sel]
        buf.decDCesaro_m2 = buf.decDCesaro_m2[:, sel_dec, sel]
        buf.decDCesaro_err = buf.decDCesaro_err[:, sel_dec, sel]

    def _fit_cesaro(self, fit=None):
        if fit is None:
            if not hasattr(self.buffer, 'fit'):
                raise Error("No fit ranges have been provided")
        else:
            fit.sort()
            self.buffer.fit = fit

    def _write_buffer(self):
        super()._write_buffer()
        self['numSample'] = self.buffer.numSample
        self['volume_err'] = self.buffer.volume_err
        self['nCorr_err'] = self.buffer.nCorr_err
        self['nDCesaro'] = self.buffer.nDCesaro
        self['nDCesaro_err'] = self.buffer.nDCesaro_err
        self['nDCesaro'].attrs['unit'] = self.buffer.nDCesaro_unit
        self['fit'] = self.buffer.fit
        self['fit'].attrs['unit'] = self.buffer.timeLags_unit
#        self['nD'] = self.buffer.nD
#        self['nD_err'] = self.buffer.nD_err
#        self['nDCesaroTotal'] = self.buffer.nDCesaroTotal
#        self['nDCesaroTotal_err'] = self.buffer.nDCesaroTotal_err
#        self['nDTotal'] = self.buffer.nDTotal
#        self['nDTotal_err'] = self.buffer.nDTotal_err

        def do_dec(dectype):
            dec_group = self.require_group(dectype.value)
            buf = getattr(self.buffer, dectype.value)
            dec_group['decCorr_err'] = buf.decCorr_err
            dec_group['decPairCount_err'] = buf.decPairCount_err

            dec_group['decDCesaro'] = buf.decDCesaro
            dec_group['decDCesaro_err'] = buf.decDCesaro_err
            dec_group['decDCesaro'].attrs['unit'] = buf.decDCesaro_unit

#            dec_group['decD'] = buf.decD
#            dec_group['decD_err'] = buf.decD_err
#            dec_group['decD'].attrs['unit'] = buf.decD_unit

        for type_ in DecType:
            if getattr(self.buffer, type_.value) is not None:
                do_dec(type_)


class Error(Exception):
    pass


def h5g_to_dict(group):
    D = {}
    for k, v in group.items():
        D[k] = v[...]
    return D


def _err_to_m2(err, n, w=None):
    """
    err: standard error of the mean
    n: number of samples
    w: average weight
    """
    if n > 1:
        if w is None:
            return np.square(err) * n * (n - 1)
        else:
            try:
                return np.square(err) * w * n * (n - 1)
            except ValueError:
                return np.square(err) * w[..., np.newaxis] * n * (n - 1)
    else:
        return np.zeros_like(err)


def _m2_to_err(m2, n, w=None):
    """
    m2: sum of squares of differences from the (current) mean
    n: number of samples
    w: average weight
    """
    if n > 1:
        if w is None:
            return np.sqrt(m2 / ((n - 1) * n))
        else:
            try:
                return np.sqrt(m2 / ((n - 1) * n * w))
            except ValueError:
                return np.sqrt(m2 / ((n - 1) * n * w[..., np.newaxis]))
    else:
        return np.full(m2.shape, np.nan)


def _get_inner_index(a, b):
    awidth = a[1] - a[0]
    bwidth = b[1] - b[0]
    if not np.isclose(awidth, bwidth):
        raise Error("Bin widths should be the same (within tolerance)")

    if (not np.isclose(np.mod(a[0] - b[0], awidth), 0.0) and
            not np.isclose(np.mod(a[0] - b[0], awidth), awidth)):
        raise Error("Bins should share a common grid.\n" +
                    "a[0]={0}, b[0]={1}, awidth={2}, (a-b)%w={3}".format(
                        a[0], b[0], awidth, np.mod(a[0] - b[0], awidth)))

    (a_begin, a_end, b_begin, b_end) = (0, len(a) - 1, 0, len(b) - 1)

    if b[0] > a[0]:
        a_begin = round((b[0] - a[0]) / awidth)
    elif b[0] < a[0]:
        b_begin = round((a[0] - b[0]) / bwidth)

    if b[-1] < a[-1]:
        a_end = round((b[-1] - a[0]) / awidth)
    elif b[-1] > a[-1]:
        b_end = round((a[-1] - b[0]) / bwidth)

    return np.s_[a_begin:a_end+1], np.s_[b_begin:b_end+1]


def new_decond(outname, samples, fit):
    with DecondFile(outname, 'x') as outfile:
        outfile._add_sample(samples, fit)
        return outfile.buffer


def extend_decond(outname, decname, samples, fit=None):
    with DecondFile(outname, 'x') as outfile:
        with DecondFile(decname) as infile:
            outfile.buffer = infile.buffer
        outfile._add_sample(samples, fit)
        return outfile.buffer


def fit_decond(outname, decname, fit):
    with DecondFile(outname, 'x') as outfile:
        with DecondFile(decname) as infile:
            outfile.buffer = infile.buffer
        outfile._fit_cesaro(fit)
        return outfile.buffer


def weighted_incremental_variance(dataweightpairs, mean=0, variance=0,
                                  numsample=0, sumweight=0):
    """
    http://www.wikiwand.com/en/Algorithms_for_calculating_variance#/Weighted_incremental_algorithm
    """
    if numsample > 1:
        m2 = variance * (numsample - 1) / numsample * sumweight
    else:
        m2 = 0

    for x, weight in dataweightpairs:
        numsample += 1
        temp = weight + sumweight
        delta = x - mean
        r = delta * weight / temp
        mean = mean + r
        m2 = m2 + sumweight * delta * r
        sumweight = temp

    variance_n = m2 / sumweight

    if numsample > 1:
        variance = variance_n * numsample / (numsample - 1)
    else:
        variance = np.full(m2.shape, np.nan)
    return mean, variance, numsample, sumweight
