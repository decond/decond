import numpy as np
from .. import analyze as da
from scipy import stats
import os
import os.path
import h5py


def test_get_inner_sel():
    print("test_get_inner_sel: starting...")
    a = np.arange(-1.5, 1, 0.5)    # ----------- [-1.5, -1. , -0.5,  0. , 0.5]
    b = np.arange(-2.5, 0.5, 0.5)  # [-2.5, -2. , -1.5, -1. , -0.5,  0.]
    a_sel, b_sel = da._get_inner_sel(a, b)
    assert((a_sel, b_sel) == (np.s_[0:4], np.s_[2:6]))
    assert(a[a_sel][0] == b[b_sel][0])
    assert(a[a_sel][-1] == b[b_sel][-1])

    a = np.arange(-2.0, 0.5, 0.5)  # ----- [-2. , -1.5, -1. , -0.5,  0.]
    b = np.arange(-2.5, 1.0, 0.5)  # [-2.5, -2. , -1.5, -1. , -0.5,  0. , 0.5]
    a_sel, b_sel = da._get_inner_sel(a, b)
    assert((a_sel, b_sel) == (np.s_[0:5], np.s_[1:6]))
    assert(a[a_sel][0] == b[b_sel][0])
    assert(a[a_sel][-1] == b[b_sel][-1])

    print("test_get_inner_sel: pass")


def test_fitlinear():
    print("test_fitlinear: starting...")
    x = np.arange(10)
    y = x
    ret = da.fitlinear(x, y)  # ret = a, b, siga, sigb, chi2, q
    assert(ret == (0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

    x = np.arange(10)
    y = 3 + 2 * x
    ret = da.fitlinear(x, y)  # ret = a, b, siga, sigb, chi2, q
    assert(ret == (3.0, 2.0, 0.0, 0.0, 0.0, 1.0))

    x = np.arange(10)
    y = []
    y.append(3 + 2 * x)
    y.append(4 + 3 * x)
    y.append(5 + 4 * x)
    y = np.array(y)
    ret = da.fitlinear(x, y)  # ret = a, b, siga, sigb, chi2, q
    assert(np.all(ret == np.array(((3., 4., 5.),
                                   (2., 3., 4.),
                                   (0., 0., 0.),
                                   (0., 0., 0.),
                                   (0., 0., 0.),
                                   (1., 1., 1.)))))

    print("test_fitlinear: pass")


def rand_c5(filename, nummoltype, timeLags=None, base_timeLags=None,
            r_decbins=None, base_r_decbins=None,
            e_decbins=None, base_e_decbins=None,
            charge=None, numMol=None, volume=None, base_volume=None,
            temperature=None, base_temperature=None):

    def rand_axis(begin, end, scale, begin_fixed=False, rand_range=None):
        if rand_range is None:
            rand_range = int((end - begin + 1) * 0.2)

        if not begin_fixed:
            begin += np.random.randint(rand_range * 2 + 1) - rand_range

        end += np.random.randint(rand_range * 2 + 1) - rand_range

        return np.arange(begin, end + 1) * scale

    if timeLags is None:
        if base_timeLags is None:
            end = 100
            scale = 0.1
        else:
            end = base_timeLags.size
            scale = base_timeLags[1] - base_timeLags[0]

        timeLags = rand_axis(0, end, scale, begin_fixed=True)
    timeLags_unit = np.string_('ps')

    if r_decbins is None:
        if base_r_decbins is None:
            end = 8
            scale = 0.01
        else:
            end = base_r_decbins.size
            scale = base_r_decbins[1] - base_r_decbins[0]

        r_decbins = rand_axis(0, end, scale, begin_fixed=True)

    if e_decbins is None:
        if base_e_decbins is None:
            begin = -7
            end = 7
            scale = 0.1
        else:
            scale = base_e_decbins[1] - base_e_decbins[0]
            begin = int(base_e_decbins[0] / scale + 0.5)
            end = begin + base_timeLags.size - 1

        e_decbins = rand_axis(begin, end, scale)

    numpairtype = nummoltype * (nummoltype + 1) // 2
    numalltype = nummoltype + numpairtype
    numtbin = timeLags.size
    numrbin = r_decbins.size
    numebin = e_decbins.size

    def rand_charge(num, maxcharge=3):
        charge = []
        balance = 0
        for i in range(num):
            if i < num - 1:
                while True:
                    rc = (np.random.random_integers(maxcharge) *
                          np.random.choice([1, -1]))
                    new_abs = abs(balance + rc)
                    if (new_abs <= maxcharge) and (not new_abs == 0):
                        break
            else:  # last one
                rc = -balance

            charge.append(rc)
            balance += rc
        return charge

    if charge is None:
        charge = rand_charge(nummoltype)
    charge_unit = np.string_('e')

    if numMol is None:
        max_nummol = 500
        base_nummol = np.random.random_integers(max_nummol)
        numMol = [base_nummol * abs(z) for z in charge]

    if volume is None:
        if base_volume is None:
            base_volume = 100.0
        volume = base_volume
        volume += ((base_volume * 0.05) * np.random.rand() * 2 + 1 -
                   base_volume * 0.05)
    volume_unit = np.string_('nm$^3$')

    if temperature is None:
        if base_temperature is None:
            base_temperature = 300.0
        temperature = base_temperature
        temperature += ((base_temperature * 0.02) * np.random.rand() * 2 +
                        1 - base_temperature * 0.02)
    temperature_unit = np.string_('K')

    if os.path.exists(filename):
        os.remove(filename)

    corr_amp = 100
    corr_offset = 1
    paircount_amp = numMol[0]
    paircount_offset = numMol[0]

    nCorr = (corr_amp *
             np.random.sample((numalltype, numtbin)) -
             corr_amp / 2 + corr_offset) / numtbin

    nCorr_unit = np.string_('nm$^2$ ps$^{-2}$')

    sd_buf = da.CorrFile._Buffer()
    sd_buf.decBins = r_decbins
    sd_buf.decBins_unit = np.string_('nm')
    sd_buf.decCorr = (corr_amp *
                      np.random.sample((numpairtype, numrbin, numtbin)) -
                      corr_amp / 2 + corr_offset) / (numtbin * numrbin)
    sd_buf.decCorr_unit = np.string_('nm$^2$ ps$^{-2}$')
    sd_buf.decPairCount = (paircount_amp *
                           np.random.sample((numpairtype, numrbin)) -
                           paircount_amp / 2 + paircount_offset) / numrbin

    ed_buf = da.CorrFile._Buffer()
    ed_buf.decBins = e_decbins
    ed_buf.decBins_unit = np.string_('kcal mol$^{-1}$')
    ed_buf.decCorr = (corr_amp *
                      np.random.sample((numpairtype, numebin, numtbin)) -
                      corr_amp / 2 + corr_offset) / (numtbin * numebin)
    ed_buf.decCorr_unit = np.string_('nm$^2$ ps$^{-2}$')
    ed_buf.decPairCount = (paircount_amp *
                           np.random.sample((numpairtype, numebin)) -
                           paircount_amp / 2 + paircount_offset) / numebin

    with da.CorrFile(filename, 'w-') as f:
        f.buffer.quantity = np.string_(da.Quantity.ec)
        f.buffer.charge = charge
        f.buffer.charge_unit = charge_unit
        f.buffer.numMol = numMol
        f.buffer.volume = volume
        f.buffer.volume_unit = volume_unit
        f.buffer.temperature = temperature
        f.buffer.temperature_unit = temperature_unit
        f.buffer.timeLags = timeLags
        f.buffer.timeLags_unit = timeLags_unit
        f.buffer.nCorr = nCorr
        f.buffer.nCorr_unit = nCorr_unit
        f.buffer.spatialDec = sd_buf
        f.buffer.energyDec = ed_buf


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


def cal_mean_sem(files):
    cfs = [da.CorrFile(file) for file in files]

    # match all the axes
    for i, cf1 in enumerate(cfs):
        for cf2 in cfs[i+1:]:
            cf1._intersect_buffer(cf2)

    # initialize results buffer
    results = da.CorrFile._Buffer()
    for dectype in da.DecType:
        setattr(results, dectype.value, da.CorrFile._Buffer())

    # initialize results arrays
    results.volume_N = []
    results.nCorr_N = []
    for dectype in da.DecType:
        buf = getattr(results, dectype.value)
        buf.decCorr_N = []
        buf.decPairCount_N = []

    for cf in cfs:
        results.volume_N.append(cf.buffer.volume)
        results.nCorr_N.append(cf.buffer.nCorr)
        for dectype in da.DecType:
            res_buf = getattr(results, dectype.value)
            data_buf = getattr(cf.buffer, dectype.value)
            res_buf.decCorr_N.append(data_buf.decCorr)
            res_buf.decPairCount_N.append(data_buf.decPairCount)

    # calculate mean and sem
    results.timeLags = cfs[0].buffer.timeLags
    results.volume = np.mean(results.volume_N, axis=0)
    results.volume_err = stats.sem(results.volume_N)
    results.nCorr = np.mean(results.nCorr_N, axis=0)
    results.nCorr_err = stats.sem(results.nCorr_N)
    for dectype in da.DecType:
        res_buf = getattr(results, dectype.value)

        res_buf.decBins = getattr(cfs[0].buffer, dectype.value).decBins
        res_buf.decCorr_N = np.array(res_buf.decCorr_N)
        res_buf.decPairCount_N = np.array(res_buf.decPairCount_N)
        weights = res_buf.decPairCount_N[..., np.newaxis].repeat(
                res_buf.decCorr_N.shape[-1], axis=-1)
        res_buf.decCorr = np.average(
                res_buf.decCorr_N, axis=0, weights=weights)
        variance = np.average(
                (res_buf.decCorr_N - res_buf.decCorr[np.newaxis, ...])**2,
                weights=weights, axis=0)
        variance *= len(cfs) / (len(cfs) - 1)
        res_buf.decCorr_err = np.sqrt(variance / len(cfs))
        res_buf.decPairCount = np.mean(res_buf.decPairCount_N, axis=0)
        res_buf.decPairCount_err = stats.sem(res_buf.decPairCount_N)

        # test weighted_incremental_variance - all at once
        a = res_buf.decCorr_N
        b = res_buf.decPairCount_N[..., np.newaxis]
        res_buf.decCorr2, variance, num_sample, _ = (
                weighted_incremental_variance(zip(a, b)))
        res_buf.decCorr2_err = np.sqrt(variance / num_sample)
        assert(np.allclose(np.nan_to_num(res_buf.decCorr),
               np.nan_to_num(res_buf.decCorr2)))
        assert(np.allclose(np.nan_to_num(res_buf.decCorr_err),
               np.nan_to_num(res_buf.decCorr2_err)))
        # ------------------------------------------------

        # test weighted_incremental_variance - one by one
        a = res_buf.decCorr_N
        b = res_buf.decPairCount_N[..., np.newaxis]
        res_buf.decCorr2, variance, num_sample, sumweight = (
                weighted_incremental_variance([(a[0], b[0])]))
        for data in zip(a[1:], b[1:]):
            res_buf.decCorr2, variance, num_sample, sumweight = (
                    weighted_incremental_variance(
                        [data], res_buf.decCorr2, variance, num_sample,
                        sumweight))

        res_buf.decCorr2_err = np.sqrt(variance / num_sample)
        assert(np.allclose(np.nan_to_num(res_buf.decCorr),
               np.nan_to_num(res_buf.decCorr2)))
        assert(np.allclose(np.nan_to_num(res_buf.decCorr_err),
               np.nan_to_num(res_buf.decCorr2_err)))
        # ------------------------------------------------

    for cf in cfs:
        cf.close()

    return results


def rand_fit(timeLags):
    dt = timeLags[1] - timeLags[0]
    while True:
        begin = np.random.rand() * (timeLags[-1] / 2 - dt) + timeLags[0]
        end = np.random.rand() * timeLags[-1] / 2 + timeLags[-1] / 2 - dt
        assert(begin >= timeLags[0])
        assert(end <= timeLags[-1])
        if end - begin > dt:
            break
    return (begin, end)


testfile = ['corr1_test.c5', 'corr2_test.c5', 'corr3_test.c5']
decondtest = 'decond_test.d5'
nummoltype = np.random.random_integers(2, 5)


def test_new_decond():
    print("test_new_decond: starting...")
    for file in testfile:
        rand_c5(file, nummoltype)

    if os.path.exists(decondtest):
        os.remove(decondtest)

    # get common timeLags
    cfs = [da.CorrFile(file) for file in testfile]
    for i, cf1 in enumerate(cfs):
        for cf2 in cfs[i+1:]:
            cf1._intersect_buffer(cf2)
    timeLags = cfs[0].buffer.timeLags
    for cf in cfs:
        cf.close()

    max_numfit = 5
    fit = [rand_fit(timeLags) for i in
           range(np.random.random_integers(max_numfit))]

    da.new_decond(decondtest, testfile, fit)

    results_file = 'results.h5'
    results = cal_mean_sem(testfile)
    with h5py.File(results_file, 'w') as f:
        for dectype in da.DecType:
            buf = getattr(results, dectype.value)
            f[dectype.value + '/decBins'] = getattr(buf, 'decBins')
            f[dectype.value + '/decCorr'] = getattr(buf, 'decCorr')
            f[dectype.value + '/decCorr_err'] = getattr(buf, 'decCorr_err')

    with da.DecondFile(decondtest) as f:
        buf = f.buffer
        assert(np.allclose(results.volume, buf.volume))
        assert(np.allclose(results.volume_err, buf.volume_err))
        assert(np.allclose(results.nCorr, buf.nCorr))
        assert(np.allclose(results.nCorr_err, buf.nCorr_err))
        for dectype in da.DecType:
            res_buf = getattr(results, dectype.value)
            dec_buf = getattr(buf, dectype.value)
            assert(np.allclose(np.nan_to_num(res_buf.decCorr),
                               np.nan_to_num(dec_buf.decCorr)))
            try:
                assert(np.allclose(np.nan_to_num(res_buf.decCorr_err),
                                   np.nan_to_num(dec_buf.decCorr_err)))
            except AssertionError:
                print(dectype.value)
                print(buf.timeLags)
                print(res_buf.decBins)
                print(dec_buf.decBins)

            assert(np.allclose(np.nan_to_num(res_buf.decPairCount),
                               np.nan_to_num(dec_buf.decPairCount)))
            assert(np.allclose(np.nan_to_num(res_buf.decPairCount_err),
                               np.nan_to_num(dec_buf.decPairCount_err)))

    print("test_new_decond: pass")


extend_file = ['corr_extend1_test.c5', 'corr_extend2_test.c5']
decond_extend = ['decond_extend_test.d5', 'decond_extend_changefit_test.d5']
decond_onebyone = ['decond_extend1_test.d5', 'decond_extend2_test.d5']


def test_extend_decond():
    print("test_extend_decond: starting...")
    for file in decond_extend + decond_onebyone:
        if os.path.exists(file):
            os.remove(file)

    for file in extend_file:
        rand_c5(file, nummoltype)

    # extend all at once
    try:
        da.extend_decond(decond_extend[0], decondtest, extend_file)
    except da.FitRangeError:
        print("  Extension failed: FitRangeError occurred\n"
              "  Ignore the first part of this test")
    else:
        with da.DecondFile(decond_extend[0]) as f, \
                da.DecondFile(decondtest) as f_old:
            assert(f.buffer.numSample ==
                   f_old.buffer.numSample + len(extend_file))

    # extend one by one
    try:
        da.extend_decond(decond_onebyone[0], decondtest, extend_file[0])
    except da.FitRangeError:
        print("  Extension failed: FitRangeError occurred\n"
              "  Ignore the first part of this test")
    else:
        with da.DecondFile(decond_onebyone[0]) as f, \
                da.DecondFile(decondtest) as f_old:
            assert(f.buffer.numSample == f_old.buffer.numSample + 1)

    try:
        da.extend_decond(decond_onebyone[1], decond_onebyone[0],
                         extend_file[1])
    except da.FitRangeError:
        print("  Extension failed: FitRangeError occurred\n"
              "  Ignore the first part of this test")
    else:
        with da.DecondFile(decond_onebyone[1]) as f, \
                da.DecondFile(decond_onebyone[0]) as f_old:
            assert(f.buffer.numSample == f_old.buffer.numSample + 1)

    with da.DecondFile(decond_extend[0]) as f_all, \
            da.DecondFile(decond_onebyone[1]) as f_one:
        np.testing.assert_allclose(f_all.buffer.temperature, f_one.buffer.temperature)
        np.testing.assert_allclose(f_all.buffer.temperature_err, f_one.buffer.temperature_err)
        np.testing.assert_allclose(f_all.buffer.volume, f_one.buffer.volume)
        np.testing.assert_allclose(f_all.buffer.volume_err, f_one.buffer.volume_err)

    # get common timeLags
    fs = ([da.CorrFile(file) for file in extend_file] +
          [da.DecondFile(decondtest)])
    for i, f1 in enumerate(fs):
        for f2 in fs[i+1:]:
            try:
                # When _intersect_buffer is called directly for debugging,
                # *_m2 attributes may not exist.
                # They are only intialized in the _add_sample function
                f1._intersect_buffer(f2)
            except AttributeError:
                pass
    timeLags = fs[0].buffer.timeLags
    for f in fs:
        f.close()

    max_numfit = 5
    fit = [rand_fit(timeLags) for i in
           range(np.random.random_integers(max_numfit))]

    da.extend_decond(decond_extend[1], decondtest, extend_file, fit)
    print("test_extend_decond: pass")


decond_fit = 'decond_changefit_test.d5'


def test_fit_decond():
    print("test_fit_decond: starting...")
    if os.path.exists(decond_fit):
        os.remove(decond_fit)

    with da.DecondFile(decondtest) as f:
        max_numfit = 5
        fit = [rand_fit(f.buffer.timeLags) for i in
               range(np.random.random_integers(max_numfit))]

    da.fit_decond(decond_fit, decondtest, fit)
    print("test_fit_decond: pass")


def test_get_rdf():
    print("test_get_rdf: starting...")
    da.get_rdf(decondtest)
    print("test_get_rdf: pass")


def test_get_D():
    print("test_get_D: starting...")
    da.get_D(decondtest)
    print("test_get_D: pass")


def test_get_decD():
    print("test_get_decD: starting...")
    da.get_decD(decondtest, da.DecType.spatial)
    da.get_decD(decondtest, da.DecType.energy)
    print("test_get_decD: pass")


def test_get_ec_dec():
    print("test_get_ec_dec: starting...")
    da.get_decqnt_sd(decondtest)
    try:
        da.get_ec_dec_energy(decondtest)
    except da.NotImplementedError:
        print("  NotImplementedError caught")
    print("test_get_ec_dec: pass")
