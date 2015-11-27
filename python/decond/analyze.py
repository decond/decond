import h5py
import numpy as np
import scipy.integrate as integrate
import scipy.constants as const
from scipy.special import gammainc
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
        self.buffer.temperature = self['temperature'][...]
        self.buffer.temperature_unit = self['temperature'].attrs['unit']
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

    @property
    def num_moltype(self):
        return _numtype(self.buffer.numMol)[0]

    @property
    def num_pairtype(self):
        return _numtype(self.buffer.numMol)[1]

    @property
    def num_alltype(self):
        return _numtype(self.buffer.numMol)[2]

    @property
    def zz(self):
        """
        Charge product for each component

        Eg. for NaCl, zz = [Na+^2, Cl-^2, Na+Na+, Na+Cl-, Cl-Cl-]
                         = [1, 1, 1, -1, 1]
        """
        return _zz(self.buffer.charge, self.buffer.numMol)

    @property
    def ww(self):
        """
        Weight for each component considering the double cross terms

        Eg. for NaCl, ww = [1, 1, 1, 2, 1]
        """
        ww = np.ones(self.num_alltype, dtype=np.int)
        for i in range(self.num_moltype):
            for j in range(i, self.num_moltype):
                if i == j:
                    ww[self.num_moltype +
                       _pairtype_index(i, j, self.num_moltype)] = 1
                else:
                    ww[self.num_moltype +
                       _pairtype_index(i, j, self.num_moltype)] = 2
        return ww

    @property
    def rdf(self):
        buf = self.buffer
        sbuf = getattr(buf, DecType.spatial.value)
        if buf is not None:
            return _paircount_to_rdf(
                    sbuf.decPairCount, sbuf.decBins, buf.numMol, buf.volume)
        else:
            raise Error("No spatialDec is found, so no rdf")

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

        self.buffer.nDTotalCesaro = np.sum(
                self.buffer.nDCesaro *
                (self.zz * self.ww)[:, np.newaxis], axis=0)

        self.buffer.nDTotalCesaro_unit = np.string_(
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
        self['temperature'] = self.buffer.temperature
        self['temperature'].attrs['unit'] = self.buffer.temperature_unit
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
        s_sel, n_sel = _get_inner_sel(
                self.buffer.timeLags, new_file.buffer.timeLags)

        self._shrink_corr_buffer(s_sel)
        new_file._shrink_corr_buffer(n_sel)
#        assert(np.allclose(self.buffer.timeLags, new_file.buffer.timeLags))

        def do_dec(dectype):
            sb_dec = getattr(self.buffer, dectype.value)
            nb_dec = getattr(new_file.buffer, dectype.value)
            s_sel_dec, n_sel_dec = _get_inner_sel(
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

    @property
    def fit_sel(self):
        return _fit_to_sel(self.buffer.fit, self.buffer.timeLags)

    def _read_decond_buffer(self):
        self.buffer.numSample = self['numSample'][...]
        self.buffer.volume_err = self['volume_err'][...]
        self.buffer.temperature_err = self['temperature_err'][...]
        self.buffer.nCorr_err = self['nCorr_err'][...]
        self.buffer.nDCesaro = self['nDCesaro'][...]
        self.buffer.nDCesaro_err = self['nDCesaro_err'][...]
        self.buffer.nDCesaro_unit = self['nDCesaro'].attrs['unit']
        self.buffer.nDTotalCesaro = self['nDTotalCesaro'][...]
        self.buffer.nDTotalCesaro_err = self['nDTotalCesaro_err'][...]
        self.buffer.nDTotalCesaro_unit = self['nDTotalCesaro'].attrs['unit']
        self.buffer.fit = self['fit'][...]
        self.buffer.fit_unit = self['fit'].attrs['unit']
        self.buffer.nD = self['nD'][...]
        self.buffer.nD_err = self['nD_err'][...]
        self.buffer.nD_unit = self['nD'].attrs['unit']
        self.buffer.nDTotal = self['nDTotal'][...]
        self.buffer.nDTotal_err = self['nDTotal_err'][...]
        self.buffer.nDTotal_unit = self['nDTotal'].attrs['unit']

        def do_dec(dectype):
            dec_group = self[dectype.value]
            buf = getattr(self.buffer, dectype.value)
            buf.decCorr_err = dec_group['decCorr_err'][...]
            buf.decPairCount_err = dec_group['decPairCount_err'][...]
            buf.decDCesaro = dec_group['decDCesaro'][...]
            buf.decDCesaro_err = dec_group['decDCesaro_err'][...]
            buf.decDCesaro_unit = dec_group['decDCesaro'].attrs['unit']
            buf.decD = dec_group['decD'][...]
            buf.decD_err = dec_group['decD_err'][...]
            buf.decD_unit = dec_group['decD'].attrs['unit']

        for type_ in DecType:
            if type_.value in self:
                do_dec(type_)

    def _add_sample(self, samples, fit, report):
        if not isinstance(samples, list):
            samples = [samples]

        if self.buffer.numSample == 0:
            if (report):
                print("Reading {0} of {1} corr files: {2}".format(
                    1, len(samples), samples[0]))
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
            init_Err('temperature')
            init_Err('nCorr')
            init_Err('nDCesaro')
            init_Err('nDTotalCesaro')

            def init_decErr(buf):
                num_sample = self.buffer.numSample

                buf.decCorr_m2 = np.zeros_like(buf.decCorr)
                buf.decCorr_err = _m2_to_err(
                        buf.decCorr_m2, num_sample)  # nan

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
            self.buffer.temperature_m2 = _err_to_m2(
                    self.buffer.temperature_err, self.buffer.numSample)
            self.buffer.nCorr_m2 = _err_to_m2(self.buffer.nCorr_err,
                                              self.buffer.numSample)
            self.buffer.nDCesaro_m2 = _err_to_m2(self.buffer.nDCesaro_err,
                                                 self.buffer.numSample)
            self.buffer.nDTotalCesaro_m2 = _err_to_m2(
                    self.buffer.nDTotalCesaro_err, self.buffer.numSample)

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
            with np.errstate(invalid='ignore'):
                r = delta * new_weight[..., np.newaxis] / temp[..., np.newaxis]
            mean += r
            m2 += sum_weight[..., np.newaxis] * delta * r
            sum_weight = temp

            if np.isscalar(mean):
                setattr(buf, data_name, mean)
                setattr(buf, data_name + '_m2', m2)

            setattr(buf, data_name + '_err', _m2_to_err(
                m2, self.buffer.numSample,
                sum_weight / self.buffer.numSample))

        def add_dec_data(dectype, new_buf):
            buf = getattr(new_buf, dectype.value)
            add_weighted_data('decCorr', 'decPairCount',
                              buf.decCorr, buf.decPairCount, dectype)
            add_weighted_data('decDCesaro', 'decPairCount',
                              buf.decDCesaro, buf.decPairCount, dectype)

            # Note that decPairCount must be updated last
            add_data('decPairCount', buf.decPairCount, dectype)

        for i, sample in enumerate(samples[begin:]):
            if (report):
                print("Reading {0} of {1} corr files: {2}".format(
                    i+begin+1, len(samples), sample))

            with CorrFile(sample) as f:
                self.buffer.numSample += 1
                self._intersect_buffer(f)
                f._cal_cesaro()

                add_data('volume', f.buffer.volume)
                add_data('temperature', f.buffer.temperature)
                add_data('nCorr', f.buffer.nCorr)
                add_data('nDCesaro', f.buffer.nDCesaro)
                add_data('nDTotalCesaro', f.buffer.nDTotalCesaro)

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
        self.buffer.nDTotalCesaro = self.buffer.nDTotalCesaro[..., sel]
        self.buffer.nDTotalCesaro_m2 = self.buffer.nDTotalCesaro_m2[..., sel]
        self.buffer.nDTotalCesaro_err = self.buffer.nDTotalCesaro_err[...,
                                                                      sel]

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
        buf = self.buffer

        if fit is None:
            try:
                buf.fit = np.asarray(buf.fit)
            except AttributeError:
                raise Error("No fit ranges have been provided")
            fit = buf.fit
        else:
            fit = sorted(fit, key=lambda x: x[0])
            fit = np.asarray(fit)
            buf.fit = fit

        fit_sel = self.fit_sel

        def fit_data(data_name, unit_ref_name, dectype=None):
            if dectype is not None:
                buf = getattr(self.buffer, dectype.value)
            else:
                buf = self.buffer
            num_sample = self.buffer.numSample
            timeLags = self.buffer.timeLags

            data_cesaro = getattr(buf, data_name + 'Cesaro')
            data_cesaro_err = getattr(buf, data_name + 'Cesaro_err')
            data_cesaro_std = _err_to_std(data_cesaro_err, num_sample)
            data_fit = np.empty((len(fit_sel),) + data_cesaro.shape[:-1])
            data_std = np.empty((len(fit_sel),) + data_cesaro.shape[:-1])

            if num_sample > 1:
                for i, sel in enumerate(fit_sel):
                    _, data_fit[i], _, data_std[i], _, _ = fitlinear(
                            timeLags[sel], data_cesaro[..., sel],
                            data_cesaro_std[..., sel])
            else:
                for i, sel in enumerate(fit_sel):
                    _, data_fit[i], _, data_std[i], _, _ = fitlinear(
                            timeLags[sel], data_cesaro[..., sel])

            data_err = _std_to_err(data_std, num_sample)

            setattr(buf, data_name, data_fit)
            setattr(buf, data_name + '_err', data_err)

            unit_ref = getattr(buf, unit_ref_name)
            unit_list = unit_ref.decode().split()
            unit_L2 = unit_list[0]
            unit_T = unit_list[1].split(sep='$')[0]
            unit_T_1 = "${0}^{{-1}}$".format(unit_T)
            unit = unit_L2 + ' ' + unit_T_1
            setattr(buf, data_name + '_unit', np.string_(unit))

        fit_data('nD', 'nCorr_unit')
        fit_data('nDTotal', 'nCorr_unit')

        for dectype in DecType:
            if getattr(buf, dectype.value) is not None:
                fit_data('decD', 'decCorr_unit', dectype)

    def _change_window(self, window):
        buf = self.buffer
        def window_data(dectype):
            decbuf = getattr(buf, dectype.value)
            binw = window[dectype]

            # removes trailing undividable bins
            decbuf.decBins = decbuf.decBins[:decbuf.decBins.size // binw * binw]
            decbuf.decBins = np.concatenate(np.mean(np.split(
                decbuf.decBins, decbuf.decBins.size // binw, axis=0), axis=1, keepdims=True), axis=0)
            initsize = decbuf.decBins.size * binw

            # decPairCount: [type, bins]
            decbuf.decPairCount = decbuf.decPairCount[:, :initsize]  # [type, decBins*binw]

            # decD: [fit, type, bins]  # L^2 T^-1
            decbuf.decD = decbuf.decD[:, :, :initsize]  # [fit, type, decBins*binw]
            decbuf.decD *= decbuf.decPairCount[np.newaxis, :, :]  # [fit, type, decBins*binw]
            decbuf.decD = np.array(np.split(decbuf.decD, decbuf.decBins.size, axis=2))  # [decBins, fit, type, binw]
            decbuf.decD = np.rollaxis(np.sum(decbuf.decD, axis=3), 0, 3)  # [fit, type, decBins]

            # decCorr: [type, bins, time]
            decbuf.decCorr = decbuf.decCorr[:, :initsize, :]  # [type, decBins*binw, time]
            decbuf.decCorr *= decbuf.decPairCount[:, :, np.newaxis]  # [type, decBins*binw, time]
            decbuf.decCorr = np.array(np.split(decbuf.decCorr, decbuf.decBins.size, axis=1))  # [decBins, type, binw, time]
            decbuf.decCorr = np.rollaxis(np.sum(decbuf.decCorr, axis=2), 0, 2)  # [type, decBins, time]

            # decDCesaro: [type, bins, time]
            decbuf.decDCesaro = decbuf.decDCesaro[:, :initsize, :]  # [type, decBins*binw, time]
            decbuf.decDCesaro *= decbuf.decPairCount[:, :, np.newaxis]  # [type, decBins*binw, time]
            decbuf.decDCesaro = np.array(np.split(decbuf.decDCesaro, decbuf.decBins.size, axis=1))  # [decBins, type, binw, time]
            decbuf.decDCesaro = np.rollaxis(np.sum(decbuf.decDCesaro, axis=2), 0, 2)  # [type, decBins, time]

            # window decPairCount
            decbuf.decPairCount = np.array(np.split(decbuf.decPairCount, decbuf.decBins.size, axis=1))  # [decBins, type, binw]
            decbuf.decPairCount = np.sum(decbuf.decPairCount, axis=2).T  # [type, decBins]

            # normalize again
            decbuf.decD /= decbuf.decPairCount[np.newaxis, :, :]
            decbuf.decCorr /= decbuf.decPairCount[:, :, np.newaxis]
            decbuf.decDCesaro /= decbuf.decPairCount[:, :, np.newaxis]

            # TODO: see how to decide error when binw > 1
            decbuf.decD_err = np.full(decbuf.decD.shape, np.nan)
            decbuf.decCorr_err = np.full(decbuf.decCorr.shape, np.nan)
            decbuf.decDCesaro_err = np.full(decbuf.decDCesaro.shape, np.nan)
            decbuf.decPairCount_err = np.full(decbuf.decPairCount.shape, np.nan)

        for dectype in window:
            if getattr(buf, dectype.value) is None:
                raise Error("{} does not have {} data".format(
                    self.filename, dectype.value))
            else:
                if window[dectype] > 1:
                    window_data(dectype)

    def _write_buffer(self):
        super()._write_buffer()
        self['numSample'] = self.buffer.numSample
        self['volume_err'] = self.buffer.volume_err
        self['temperature_err'] = self.buffer.temperature_err
        self['nCorr_err'] = self.buffer.nCorr_err
        self['nDCesaro'] = self.buffer.nDCesaro
        self['nDCesaro_err'] = self.buffer.nDCesaro_err
        self['nDCesaro'].attrs['unit'] = self.buffer.nDCesaro_unit
        self['nDTotalCesaro'] = self.buffer.nDTotalCesaro
        self['nDTotalCesaro_err'] = self.buffer.nDTotalCesaro_err
        self['nDTotalCesaro'].attrs['unit'] = self.buffer.nDTotalCesaro_unit
        self['fit'] = self.buffer.fit
        self['fit'].attrs['unit'] = self.buffer.timeLags_unit
        self['nD'] = self.buffer.nD
        self['nD_err'] = self.buffer.nD_err
        self['nD'].attrs['unit'] = self.buffer.nD_unit
        self['nDTotal'] = self.buffer.nDTotal
        self['nDTotal_err'] = self.buffer.nDTotal_err
        self['nDTotal'].attrs['unit'] = self.buffer.nDTotal_unit

        def do_dec(dectype):
            dec_group = self.require_group(dectype.value)
            buf = getattr(self.buffer, dectype.value)
            dec_group['decCorr_err'] = buf.decCorr_err
            dec_group['decPairCount_err'] = buf.decPairCount_err

            dec_group['decDCesaro'] = buf.decDCesaro
            dec_group['decDCesaro_err'] = buf.decDCesaro_err
            dec_group['decDCesaro'].attrs['unit'] = buf.decDCesaro_unit

            dec_group['decD'] = buf.decD
            dec_group['decD_err'] = buf.decD_err
            dec_group['decD'].attrs['unit'] = buf.decD_unit

        for type_ in DecType:
            if getattr(self.buffer, type_.value) is not None:
                do_dec(type_)


class Error(Exception):
    pass


class FitRangeError(Error):
    pass


class NotImplementedError(Error):
    pass


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
            # it is important to explicitly extend the dimension of weight, w,
            # by [..., np.newaxis], otherwise, numpy will broadcast it
            # wrongly when the dimension of m2 and w are [N, N, N] and [N, N],
            # respectivly.
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
            # it is important to explicitly extend the dimension of weight, w,
            # by [..., np.newaxis], otherwise, numpy will broadcast it
            # wrongly when the dimension of m2 and w are [N, N, N] and [N, N],
            # respectivly.
            return np.sqrt(m2 / ((n - 1) * n * w[..., np.newaxis]))
    else:
        return np.full(m2.shape, np.nan)


def _err_to_var(err, n):
    """
    err: standard error of the mean
    n: number of samples
    """
    if n > 1:
        err = np.square(err) * n
        err[err == 0] = np.nan
        return err
    else:
        return np.full(err.shape, np.nan)


def _err_to_std(err, n):
    """
    err: standard error of the mean
    n: number of samples
    """
    if n > 1:
        return err * np.sqrt(n)
    else:
        return np.full(err.shape, np.nan)


def _std_to_err(std, n):
    """
    std: standard deviation
    n: number of samples
    """
    if n > 1:
        return std / np.sqrt(n)
    else:
        return np.full(std.shape, np.nan)


def _get_inner_sel(a, b):
    """
    Return the intersection of a and b in terms of np.s_ with respect to
    a and b, respectively.
    """
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


def _pairtype_index(moltype1, moltype2, num_moltype):
    """
    Return pairtype from two moltypes
              c
        | 0  1  2  3
      --+------------
      0 | 0  1  2  3
        |
      1 |    4  5  6
    r   |
      2 |       7  8
        |
      3 |          9

      index(r, c) = r * n + c - r * (r + 1) / 2
      where n = size(c) = size(r), r <= c
    """
    r = min(moltype1, moltype2)
    c = max(moltype1, moltype2)
    return r * num_moltype + c - r * (r + 1) / 2


def _fit_to_sel(fit, timelags):
    """
    Return the selection of type np.s_ corresponding to the fit range
    """
    sel = []
    fit = np.asarray(fit)
    if fit.ndim != 2 or fit.shape[1] != 2:
        raise Error("fit should be of shape (N, 2)")

    if np.any(fit < timelags[0]) or np.any(fit > timelags[-1]):
        raise FitRangeError(
                "fit range is out of timelags range...\n" +
                "timelags[0]:{0}, timelags[-1]:{1}\n".format(
                    timelags[0], timelags[-1]) + "fit:{0}".format(fit))

    dt = timelags[1] - timelags[0]
    for fit_ in fit:
        if fit_[1] <= fit_[0]:
            raise Error("Unreasonable fit, end <= begin: {0}".format(fit_))
        begin, end = (fit_ / dt).astype(int)
        sel.append(np.s_[begin:end])
    return sel


def fitlinear(x, y, sig=None):
    """
    Given a set of data points x, y with individual standard deviations sig,
    fit them to a straight line y = a + bx by minimizing chi2.

    Returned are a,b and their respective probable uncertainties siga and sigb,
    the chi-square chi2, and the goodness-of-fit probability q (that the fit
    would have chi2 this large or larger).
    If no sig is provided, the standard deviations are assumed to be
    unavailable: q is returned as 1.0 and the normalization of chi2
    is to unit standard deviation on all points.

    Parameters
    x: 1-dimension
    y: 1 or more dimensions, fit only the last axis
    sig (optional): same dimension as y

    Return:
    a, b, siga, sigb, chi2, q
    """
    if x.ndim != 1:
        raise Error("x must be one dimensional. x.ndim={0}".format(x.ndim))
    if x.size != y.shape[-1]:
        raise Error("lengths of the last dimension of x and y do not match\n"
                    "x.size={0}, y.shape[-1]={1}".format(x.size, y.shape[-1]))
    if sig is not None:
        sig2 = sig**2
        sig2[sig2 == 0] = np.nan  # TODO: not sure if it is a good solution
        wt = 1 / sig2  # y.shape
        ss = np.sum(wt, axis=-1)  # y.shape[:-1]
        sx = np.sum(x * wt, axis=-1)  # y.shape[:-1]
        sy = np.sum(y * wt, axis=-1)  # y.shape[:-1]
    else:
        sx = np.sum(x, axis=-1)  # scalar
        sy = np.sum(y, axis=-1)  # y.shape[:-1]
        ss = x.size  # scalar
    sxoss = sx / ss  # y.shape[:-1] or scalar
    if sig is not None:
        t = (x - sxoss[..., np.newaxis]) / sig  # y.shape
        st2 = np.sum(t * t, axis=-1)  # y.shape[:-1]
        b = np.sum(t * y / sig, axis=-1)  # y.shape[:-1]
    else:
        t = x - sxoss  # x.shape
        st2 = np.sum(t * t, axis=-1)  # scalar
        b = np.sum(t * y, axis=-1)  # y.shape[:-1]

    b = b / st2  # y.shape[:-1]
    a = (sy - sx * b) / ss  # y.shape[:-1]
    siga = np.sqrt((1. + sx * sx / (ss * st2)) / ss)  # y.shape[:-1] or scalar
    sigb = np.sqrt(1. / st2)  # y.shape[:-1] or scalar
    q = np.ones(y.shape[:-1])
    if sig is None:
        chi2 = np.sum((y - a[..., np.newaxis] - b[..., np.newaxis] * x)**2,
                      axis=-1)  # y.shape[:-1]
        sigdat = np.sqrt(chi2 / (x.size - 2))  # y.shape[:-1]
        siga = siga * sigdat  # y.shape[:-1]
        sigb = sigb * sigdat  # y.shape[:-1]
    else:
        chi2 = np.sum(
                ((y - a[..., np.newaxis] - b[..., np.newaxis] * x) / sig)**2,
                axis=-1)  # y.shape[:-1]
        if x.size > 2:
            q = gammainc(0.5 * (x.size - 2), 0.5 * chi2)  # y.shape[:-1]

    return a, b, siga, sigb, chi2, q


def _nummolpair(nummol):
    return np.array([n1 * (n2-1) if e1 == e2 else n1*n2
                     for (e1, n1) in enumerate(nummol)
                     for (e2, n2) in enumerate(nummol) if e2 >= e1])


def _paircount_to_rdf(paircount, rbins, nummol, volume):
    v2 = volume**2
    pair_density = _nummolpair(nummol) / v2

    l_half = volume**(1/3) / 2
    dr = rbins[1] - rbins[0]
    dv = 4 * np.pi * dr * rbins**2
    dv[dv == 0] = np.nan
    dvsim = dv.copy()
    filter = (rbins > l_half) & (rbins < np.sqrt(2) * l_half)
    dvsim[filter] = 4 * np.pi * dr * rbins[filter] * (3 * l_half -
                                                      2 * rbins[filter])
    filter = (rbins >= np.sqrt(2) * l_half)
    dvsim[filter] = np.nan
    rho_vdvsim = paircount
    rho_v = rho_vdvsim / dvsim
    # rho_dvsim = rho_vdvsim / volume
    rho = rho_v / volume

    g = rho / pair_density[:, np.newaxis]
    return g


def _paircount_to_edf(paircount, ebins, volume):
    de = ebins[1] - ebins[0]
    rho_vde = paircount
    rho_de = rho_vde / volume
    edf = rho_de / de
    # normalize
    # edf = rho_de / integrate.trapz(rho_de)[..., np.newaxis]
    return edf  # L^3 E^-1


def _numtype(nummol):
    """
    Return: num_moltype, num_pairtype, num_alltype
    """
    num_moltype = nummol.size
    num_pairtype = num_moltype * (num_moltype + 1) // 2
    num_alltype = num_moltype + num_pairtype
    return num_moltype, num_pairtype, num_alltype


def _zz(charge, nummol):
    num_moltype, _, num_alltype = _numtype(nummol)
    zz = np.empty(num_alltype, dtype=np.int)
    for i in range(num_moltype):
        zz[i] = charge[i]**2
        for j in range(i, num_moltype):
            zz[num_moltype + _pairtype_index(i, j, num_moltype)] = (
                    charge[i] * charge[j])
    return zz


def _nD_to_ec_const(temperature, volume):
    beta = 1 / (const.k * temperature)
    nD_to_ec = beta / volume
    nD_to_ec *= const.e**2 / const.nano**3 * (const.nano**2 / const.pico)
    return nD_to_ec


def get_fit(decname):
    """
    Return fit, fit_unit
    """
    with h5py.File(decname, 'r') as f:
        fit = f['fit'][...]
        fit *= const.pico
        fit_unit = "s"
    return fit, fit_unit


def get_rdf(decname):
    """
    Return rdf, rbins, rbins_unit
    """
    with h5py.File(decname, 'r') as f:
        gid = f[DecType.spatial.value]
        rbins = gid['decBins'][...]
        rdf = _paircount_to_rdf(gid['decPairCount'][...], rbins,
                                f['numMol'][...], f['volume'][...])
        rbins *= const.nano
        rbins_unit = "m"

        return rdf, rbins, rbins_unit


def get_edf(decname):
    """
    Return edf, ebins, ebins_unit
    """
    with h5py.File(decname, 'r') as f:
        gid = f[DecType.energy.value]
        ebins = gid['decBins'][...] * const.calorie  # kJ mol^-1
        volume = f['volume'][...] * const.nano**3  # m^3
        edf = _paircount_to_edf(gid['decPairCount'][...], ebins, volume)
        edf_unit = "m$^{-3}$ kJ$^{-1} mol"
        ebins_unit = "kJ mol$^{-1}$"
        return edf, edf_unit, ebins, ebins_unit


def get_diffusion(decname):
    """
    Return diffusion, diffusion_err, diffusion_unit, fit, fit_unit
    """
    with h5py.File(decname, 'r') as f:
        nummol = f['numMol'][...]
        num_moltype, _, _ = _numtype(nummol)
        nD = f['nD'][:, :num_moltype]
        nD_err = f['nD_err'][:, :num_moltype]

        diffusion = nD / nummol  # L^2 T^-1  [fit, num_moltype]
        diffusion_err = nD_err / nummol
        cc = const.nano**2 / const.pico
        diffusion *= cc
        diffusion_err *= cc
        diffusion_unit = "m$^2$ s$^{-1}$"

    fit, fit_unit = get_fit(decname)
    return diffusion, diffusion_err, diffusion_unit, fit, fit_unit


def get_decD(decname, dectype):
    """
    Return decD, decD_err, decD_unit, decBins, decBins_unit, fit, fit_unit
    """
    with h5py.File(decname, 'r') as f:
        gid = f[dectype.value]
        decD = gid['decD'][...]  # L^2 T^-1
        decD_err = gid['decD_err'][...]
        decBins = gid['decBins'][...]
        ccD = const.nano**2 / const.pico
        decD *= ccD
        decD_err *= ccD
        decD_unit = "m$^2$ s$^{-1}$"
        if dectype is DecType.spatial:
            decBins *= const.nano
            decBins_unit = "m"
        elif dectype is DecType.energy:
            decBins *= const.calorie
            decBins_unit = "kJ mol$^{-1}$"

    fit, fit_unit = get_fit(decname)
    return decD, decD_err, decD_unit, decBins, decBins_unit, fit, fit_unit


def get_ec(decname):
    """
    Return
    ec_total, ec_total_err,
    ec[:num_moltype], ec_err[:num_moltype], ec_unit,
    fit, fit_unit
    """
    with h5py.File(decname, 'r') as f:
        volume = f['volume'][...]
        temperature = f['temperature'][...]
        nDTotal = f['nDTotal'][...]
        nDTotal_err = f['nDTotal_err'][...]
        nD = f['nD'][...]
        nD_err = f['nD_err'][...]
        nummol = f['numMol'][...]
        charge = f['charge'][...]
        zz = _zz(charge, nummol)
        cc = _nD_to_ec_const(temperature, volume)
        ec_total = nDTotal * cc
        ec_totol_err = nDTotal_err * cc
        ec = nD * zz * cc
        ec_err = nD_err * abs(zz) * cc
        ec_unit = "S m$^{-1}$"

    fit, fit_unit = get_fit(decname)
    return ec_total, ec_totol_err, ec, ec_err, ec_unit, fit, fit_unit


def _symmetrize_array(arr, decBins, center=0, axis=-1):
    center_idx = np.where(decBins==center)[0][0]
    num_left = center_idx
    num_right = decBins.size - num_left - 1
    bw = decBins[1] - decBins[0]
    if num_left > num_right:
        # postpad zero
        numpad = num_left - num_right
        arr = np.append(arr, np.zeros(arr.shape[:-1]+(numpad,)), axis=axis)
        decBins = np.hstack((decBins, [decBins[-1] + (i+1) * bw for i in range(numpad)]))
    elif num_left < num_right:
        # prepad zero
        numpad = num_right - num_left
        arr = np.insert(arr, [0]*numpad, 0., axis=axis)
        decBins = np.hstack((list(reversed([decBins[0] - (i+1) * bw for i in range(numpad)])), decBins))

    return arr, decBins


def get_ec_dec(decname, dectype, sep_nonlocal=True, nonlocal_ref=None,
               avewidth=None):
    """
    Return ec_dec, ec_dec_unit, decBins, decBins_unit, fit, fit_unit
    """
    with h5py.File(decname, 'r') as f:
        gid = f[dectype.value]
        nummol = f['numMol'][...]
        charge = f['charge'][...]
        volume = f['volume'][...]
        temperature = f['temperature'][...]
        num_moltype, _, _ = _numtype(nummol)
        zz = _zz(charge, nummol)
        nD = f['nD'][...]
        decD = gid['decD'][...]  # L^2 T^-1
        decBins = gid['decBins'][...]
        beta = 1 / (const.k * temperature)
        paircount = gid['decPairCount'][...]
        bw = decBins[1] - decBins[0]

        if dectype is DecType.energy:
            decD, _ = _symmetrize_array(decD, decBins)
            paircount, decBins = _symmetrize_array(paircount, decBins)

        if sep_nonlocal:
            if dectype is DecType.spatial:
                if nonlocal_ref is None:
                    nonlocal_ref_idx = decBins.size / np.sqrt(3)
                else:
                    nonlocal_ref_idx = int(nonlocal_ref / bw)
                if avewidth is None:
                    avewidth = 0.25
                avewidth_idx = int(avewidth / bw)
            elif dectype is DecType.energy:
                raise NotImplementedError(
                        "ec_dec with sep_nonlocal is not implemented "
                        "for DecType.energy yet")
            else:
                raise Error("Impossible DecType...")
            decD_nonlocal = np.mean(
                    decD[:, :,
                         nonlocal_ref_idx - avewidth_idx:
                         nonlocal_ref_idx + avewidth_idx],
                    axis=-1)
            ec_nonlocal = (_nummolpair(nummol) / volume *
                           decD_nonlocal *
                           zz[num_moltype:] * beta)
        else:
            decD_nonlocal = np.zeros(decD.shape[:-1])
            ec_nonlocal = np.zeros_like(decD_nonlocal)

        ec_local = (paircount / volume *
                    (decD - decD_nonlocal[:, :, np.newaxis]) *
                    zz[num_moltype:, np.newaxis] * beta)
        ec_local[np.isnan(ec_local)] = 0

        if dectype is DecType.energy:
            ec_local = (integrate.cumtrapz(ec_local[..., :(decBins.size+1)/2], initial=0) +
                        integrate.cumtrapz(ec_local[..., -1:(decBins.size-1)/2-1:-1], initial=0))
            ec_nonlocal = np.zeros(ec_local.shape[:-1])
            decBins = decBins[:(decBins.size+1)/2]
        else:
            ec_local = integrate.cumtrapz(ec_local, initial=0)

        cc = (1 / const.nano**3 *
              const.nano**2 / const.pico *
              const.e**2)
        ec_nonlocal *= cc
        ec_local *= cc
        ec_dec_unit = "S m$^{-1}$"
        if dectype is DecType.spatial:
            decBins *= const.nano
            decBins_unit = "m"
        elif dectype is DecType.energy:
            decBins *= const.calorie
            decBins_unit = "kJ mol$^{-1}$"

        ec_auto = (nD[:, :num_moltype] * zz[:num_moltype] *
                   _nD_to_ec_const(temperature, volume))
        ec_dec = ec_auto[:, :, np.newaxis] * np.ones_like(decBins)

        # ec_dec = ec_auto + ec_cross = ec_auto + (ec_local + ec_nonlocal)
        for r in range(num_moltype):
            for c in range(r, num_moltype):
                idx = _pairtype_index(r, c, num_moltype)
                ec_dec[:, r] += (ec_local[:, idx] +
                                 ec_nonlocal[:, idx, np.newaxis])
                if (r != c):
                    ec_dec[:, c] += (ec_local[:, idx] +
                                     ec_nonlocal[:, idx, np.newaxis])

    fit, fit_unit = get_fit(decname)
    return ec_dec, ec_dec_unit, decBins, decBins_unit, fit, fit_unit


def get_normalize_paircount(decname, dectype):
    with h5py.File(decname, 'r') as f:
        gid = f[dectype.value]
        paircount = gid['decPairCount'][...]
    return paircount / integrate.trapz(paircount)[..., np.newaxis]


def get_ec_dec_energy(decname, sep_nonlocal=True, threshold=0):
    """
    Return ec_dec_cross_IL, ec_dec_cross_IL_unit, decBins, decBins_unit, fit, fit_unit
    """
    dectype = DecType.energy
    with h5py.File(decname, 'r') as f:
        gid = f[dectype.value]
        nummol = f['numMol'][...]
        charge = f['charge'][...]
        volume = f['volume'][...]
        fit = f['fit'][...]
        temperature = f['temperature'][...]
        decD = gid['decD'][...]  # L^2 T^-1
        decBins = gid['decBins'][...]
        paircount = gid['decPairCount'][...]

    beta = 1 / (const.k * temperature)
    zz = _zz(charge, nummol)
    num_moltype, num_pairtype, _ = _numtype(nummol)

    def sep_at_end():
        nonlocal_idx = np.empty(decD.shape[:-1])
        for f in range(nonlocal_idx.shape[0]):
            for t in range(nonlocal_idx.shape[1]):
                if zz[num_moltype + t] > 0:
                    nonlocal_idx[f, t] = np.nonzero(np.invert(np.isnan(decD[f, t, :])))[0][0]
                else:
                    nonlocal_idx[f, t] = np.nonzero(np.invert(np.isnan(decD[f, t, :])))[0][-1]
        for f in range(nonlocal_idx.shape[0]):
            for t in range(nonlocal_idx.shape[1]):
                decD_nonlocal[f, t] = decD[f, t, nonlocal_idx[f, t]]
        return decD_nonlocal

    def sep_at_edf_max():
        nonlocal_idx = np.empty(decD.shape[:-1])
        for f in range(nonlocal_idx.shape[0]):
            nonlocal_idx[f, :] = np.argmax(paircount, axis=-1)
        for f in range(nonlocal_idx.shape[0]):
            for t in range(nonlocal_idx.shape[1]):
                decD_nonlocal[f, t] = decD[f, t, nonlocal_idx[f, t]]
        return decD_nonlocal

    decD_nonlocal = np.empty(decD.shape[:-1])
    if sep_nonlocal:
        #decD_nonlocal = sep_at_end()
        decD_nonlocal = sep_at_edf_max()
    else:
        decD_nonlocal = np.zeros(decD.shape[:-1])

    ec_nonlocal = (
            integrate.trapz(paircount) / volume * decD_nonlocal *
            zz[np.newaxis, num_moltype:] * beta)
    ec_local = (
            paircount / volume * (decD - decD_nonlocal[:, :, np.newaxis]) *
            zz[num_moltype:, np.newaxis] * beta)

    norm_paircount = get_normalize_paircount(decname, dectype)
    norm_paircount, _ = _symmetrize_array(norm_paircount, decBins)
    ec_local, decBins = _symmetrize_array(ec_local, decBins)

    ec_local[np.isnan(ec_local)] = 0
    # filter with threshold
    # ec_local_masked = np.ma.masked_where(
            # np.ones_like(ec_local) * norm_paircount[np.newaxis, ...] < threshold,
            # ec_local)

    np.place(ec_local,
             np.ones_like(ec_local) * norm_paircount[np.newaxis, ...] < threshold,
             0)

    # reverse the integrate direction for zz > 0 component
    for i, czz in enumerate(zz[num_moltype:]):
        if czz > 0:
            ec_local[:, i, :] = ec_local[:, i, ::-1]

    ec_local = integrate.cumtrapz(ec_local, initial=0)

    ec_dec_cross_IL = ec_local + ec_nonlocal[:, :, np.newaxis]
    ec_dec_cross_I = np.zeros((fit.shape[0], num_moltype, ec_dec_cross_IL.shape[-1]))
    for r in range(num_moltype):
        for c in range(r, num_moltype):
            idx = _pairtype_index(r, c, num_moltype)
            ec_dec_cross_I[:, r] += ec_dec_cross_IL[:, idx]
            if c > r:
                ec_dec_cross_I[:, c] += ec_dec_cross_IL[:, idx]

    cc = (1 / const.nano**3 *
          const.nano**2 / const.pico *
          const.e**2)
    ec_dec_cross_IL *= cc
    ec_dec_cross_IL_unit = "S m$^{-1}$"
    ec_dec_cross_I *= cc
    ec_dec_cross_I_unit = "S m$^{-1}$"

    decBins *= const.calorie
    decBins_unit = "kJ mol$^{-1}$"

    fit, fit_unit = get_fit(decname)

    return (ec_dec_cross_I, ec_dec_cross_I_unit,
            ec_dec_cross_IL, ec_dec_cross_IL_unit,
            decBins, decBins_unit,
            fit, fit_unit)


def new_decond(outname, samples, fit, report=True):
    with DecondFile(outname, 'w-') as outfile:
        outfile._add_sample(samples, fit, report)
        return outfile.buffer


def extend_decond(outname, decname, samples, fit=None, report=True):
    with DecondFile(outname, 'w-') as outfile:
        if (report):
            print("Reading decond file: {0}".format(decname))
        with DecondFile(decname) as infile:
            outfile.buffer = infile.buffer
        outfile._add_sample(samples, fit, report)
        return outfile.buffer


def fit_decond(outname, decname, fit, report=True):
    with DecondFile(outname, 'w-') as outfile:
        if (report):
            print("Reading decond file: {0}".format(decname))
        with DecondFile(decname) as infile:
            outfile.buffer = infile.buffer
        outfile._fit_cesaro(fit)
        return outfile.buffer


def window_decond(outname, decname, window, report=True):
    with DecondFile(outname, 'w-') as outfile:
        if (report):
            print("Reading decond file: {0}".format(decname))
        with DecondFile(decname) as infile:
            outfile.buffer = infile.buffer
        outfile._change_window(window)
        return outfile.buffer


def report_decond(decname):
    print()

    diffusion, diffusion_err, diffusion_unit, fit, fit_unit = \
        get_diffusion(decname)
    print("Diffusion")
    print("=========")
    print("{:<15} {:<}".format(
        'Fit (ps)', 'Diffusion (' + diffusion_unit + ')'))
    for i in range(len(fit)):
        print("{:<15}    {:<}".format(str(fit[i] * 1e12), str(diffusion[i])))
        print("{:<15} +/-{:<}".format('', str(diffusion_err[i])))
        print()

    print()

    ec_total, ec_total_err, ec, ec_err, ec_unit, fit, fit_unit = \
        get_ec(decname)
    print("Electrical conductivity")
    print("=======================")
    print("{:<15} {:<15}".format(
        'Fit (ps)', 'Electrical conductivity (' + ec_unit + ')'))
    for i in range(len(fit)):
        print("{:<15}    {:<}".format(str(fit[i] * 1e12), str(ec[i])))
        print("{:<15} +/-{:<}".format('', str(ec_err[i])))
        print("{:<15} Total: {:<} +/- {:<}".format(
            '', str(ec_total[i]), str(ec_total_err[i])))
        print()

    print()
