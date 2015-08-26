import numpy as np
import decond.analyzer as da
import os
import os.path


def test_get_inner_index():
    a = np.arange(-1.5, 1, 0.5)    # ----------- [-1.5, -1. , -0.5,  0. , 0.5]
    b = np.arange(-2.5, 0.5, 0.5)  # [-2.5, -2. , -1.5, -1. , -0.5,  0.]
    a_sel, b_sel = da._get_inner_index(a, b)
    assert((a_sel, b_sel) == (np.s_[0:4], np.s_[2:6]))
    assert(a[a_sel][0] == b[b_sel][0])
    assert(a[a_sel][-1] == b[b_sel][-1])

    a = np.arange(-2.0, 0.5, 0.5)  # ----- [-2. , -1.5, -1. , -0.5,  0.]
    b = np.arange(-2.5, 1.0, 0.5)  # [-2.5, -2. , -1.5, -1. , -0.5,  0. , 0.5]
    a_sel, b_sel = da._get_inner_index(a, b)
    assert((a_sel, b_sel) == (np.s_[0:5], np.s_[1:6]))
    assert(a[a_sel][0] == b[b_sel][0])
    assert(a[a_sel][-1] == b[b_sel][-1])

    print("_get_inner_index: pass")

test_get_inner_index()


def gen_rand_c5(filename, nummoltype, timeLags=None, base_timeLags=None,
                r_decbins=None, base_r_decbins=None,
                e_decbins=None, base_e_decbins=None,
                charge=None, numMol=None, volume=None, base_volume=None):

    def rand_axis(begin, end, scale, begin_fixed=False, rand_range=None):
        if rand_range is None:
            rand_range = int((end - begin + 1) * 0.2)

        if not begin_fixed:
            begin += np.random.randint(rand_range * 2 + 1) - rand_range

        end += np.random.randint(rand_range * 2 + 1) - rand_range

        return np.arange(begin, end + 1) * scale

    if timeLags is None:
        if base_timeLags is None:
            end = 10
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
            end = base_timeLags.size
            scale = base_timeLags[1] - base_timeLags[0]

        r_decbins = rand_axis(0, end, scale, begin_fixed=True)

    if e_decbins is None:
        if base_e_decbins is None:
            begin = -7
            end = 7
            scale = 0.1
        else:
            end = base_timeLags.size
            scale = base_timeLags[1] - base_timeLags[0]

        e_decbins = rand_axis(begin, end, scale)

    numpairtype = nummoltype * (nummoltype + 1) / 2
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

    if os.path.exists(filename):
        os.remove(filename)

    corr_amp = 100
    corr_offset = 1
    paircount_amp = numMol[0]
    paircount_offset = numMol[0]

    nCorr = (corr_amp *
             np.random.sample((numalltype, numtbin)) -
             corr_amp / 2 + corr_offset) / numrbin

    nCorr_unit = np.string_('nm$^2$ ps$^{-2}$')

    sd_buf = da.CorrFile._Buffer()
    sd_buf.decBins = r_decbins
    sd_buf.decBins_unit = np.string_('nm')
    sd_buf.decCorr = (corr_amp *
                      np.random.sample((numpairtype, numrbin, numtbin)) -
                      corr_amp / 2 + corr_offset) / numrbin
    sd_buf.decCorr_unit = np.string_('nm$^2$ ps$^{-2}$')
    sd_buf.decPairCount = (paircount_amp *
                           np.random.sample((numpairtype, numrbin)) -
                           paircount_amp / 2 + paircount_offset) / numrbin

    ed_buf = da.CorrFile._Buffer()
    ed_buf.decBins = e_decbins
    ed_buf.decBins_unit = np.string_('kcal mol$^{-1}$')
    ed_buf.decCorr = (corr_amp *
                      np.random.sample((numpairtype, numebin, numtbin)) -
                      corr_amp / 2 + corr_offset) / numebin
    ed_buf.decCorr_unit = np.string_('nm$^2$ ps$^{-2}$')
    ed_buf.decPairCount = (paircount_amp *
                           np.random.sample((numpairtype, numebin)) -
                           paircount_amp / 2 + paircount_offset) / numebin

    with da.CorrFile(filename, 'x') as f:
        f.buffer.charge = charge
        f.buffer.charge_unit = charge_unit
        f.buffer.numMol = numMol
        f.buffer.volume = volume
        f.buffer.volume_unit = volume_unit
        f.buffer.timeLags = timeLags
        f.buffer.timeLags_unit = timeLags_unit
        f.buffer.nCorr = nCorr
        f.buffer.nCorr_unit = nCorr_unit
        f.buffer.spatialDec = sd_buf
        f.buffer.energyDec = ed_buf


testfile = ['rand1.c5', 'rand2.c5', 'rand3.c5']
for file in testfile:
    gen_rand_c5(file, 3)


def test_cal_decond():
    da.cal_decond('decond_test.d5', testfile)
    print("cal_decond: pass")

test_cal_decond()


# /charge                  Dataset {2}
# /energyDec               Group
# /energyDec/decBins       Dataset {185}
# /energyDec/decCorr       Dataset {3, 185, 11}
# /energyDec/decPairCount  Dataset {3, 185}
# /nCorr                   Dataset {5, 11}
# /numMol                  Dataset {2}
# /spatialDec              Group
# /spatialDec/decBins      Dataset {411}
# /spatialDec/decCorr      Dataset {3, 411, 11}
# /spatialDec/decPairCount Dataset {3, 411}
# /timeLags                Dataset {11}
# /volume                  Dataset {SCALAR}
