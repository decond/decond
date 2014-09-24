import matplotlib.cm as cm
windows = np.arange(0, 0.51, 0.05)[1::2]
#windows = [0.25]
for fit in sortedKeys:
  sdD[fit] /= Const.D2AA2_ps
  plt.figure() 
  plt.gca().set_color_cycle(cm.rainbow(np.linspace(0, 1, len(windows))))

  for aveWidth in windows:
    aveWidth /= dr
    # determine D_IL(\infty)
    sdDInfty = {}
    sigNonLocal = {}
    sdDInfty[fit] = np.mean(sdD[fit][:, halfCellIndex - aveWidth: halfCellIndex + aveWidth], axis=1)
  #    sdDInfty[fit] = np.zeros_like(sdDInfty[fit])
    sigNonLocal[fit] = numMolPair / (volume * Const.nm**3) * sdDInfty[fit] * Const.nm**2 / Const.ps *\
                       zzCross * Const.beta * Const.basicCharge**2

    sigIL = {}
  #    sigIL[fit] = rho_dvsim / Const.nm**3 * sdD[fit] * Const.nm**2 / Const.ps * \
  #                 zzCross[:, np.newaxis] * Const.beta * Const.basicCharge**2
    sigIL[fit] = rho_dvsim / Const.nm**3 * (sdD[fit] - sdDInfty[fit][:, np.newaxis]) *\
                 Const.nm**2 / Const.ps * zzCross[:, np.newaxis] * Const.beta * Const.basicCharge**2
    sigIL[fit][np.isnan(sigIL[fit])] = 0
    sigIL[fit] = integrate.cumtrapz(sigIL[fit], initial=0)

    sigI = {}
    sigI[fit] = sigAutoI[fit][:, np.newaxis] * np.ones_like(rBins)
    for r in range(numIonTypes):
      for c in range(r, numIonTypes):
        sigI[fit][r] += sigIL[fit][zipIndexPair2(r,c, numIonTypes)] +\
                        sigNonLocal[fit][zipIndexPair2(r,c, numIonTypes)]
        if (r != c):
          sigI[fit][c] += sigIL[fit][zipIndexPair2(r,c, numIonTypes)] +\
                          sigNonLocal[fit][zipIndexPair2(r,c, numIonTypes)]

    # plot sig
    for i, sig in enumerate([sigI[fit][0]]):
      plt.plot(rBins, sig, label=label[i] + ', win ' + str(aveWidth*dr*10) + r' $\AA$')
      plt.legend(loc='upper right')
      plt.ylabel(r"$\sigma_I(\lambda)$  (S m$^{-1}$)")
      plt.xlabel(r"$r$  ($\AA$)")

    plt.title(fit)
#    plt.xlim(xmax=cellLengthHalf*Const.nm2AA)
