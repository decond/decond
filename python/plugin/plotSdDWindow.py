def plotSdDWindow(window):
  """
  Plots D_IL(r) with windowed averages for determining D_IL(\infty)
  """
  idxPerWindow = window // (dr * Const.nm2AA)
  numWindow = rBins.size // idxPerWindow
  rBinsWindow = np.mean(np.split(rBins[:idxPerWindow * numWindow], numWindow),
                       axis=1)
  #split => [numWindow, idxPerWindow], mean => [numWindow]
  sdDWindow = {}
  for fit in sdD:
    sdDWindow[fit] = np.mean(np.split(sdD[fit][:, :idxPerWindow * numWindow],
                                     numWindow, axis=1), axis=2).T
    #split => [numWindow, type, idxPerWindow]
    #mean => [numWindow, type] 
    #transpose => [type, numWindow] 

  for fitKey in sorted(sdD, key=lambda x:x.split(sep='-')[0]):
    plt.figure()
    ax = plt.gca()
    if (args.color is not None):
      ax.set_color_cycle(args.color[numIonTypes:])
    for i, D in enumerate(sdDWindow[fitKey]):
      plt.plot(rBinsWindow, D, label=label[numIonTypes + i],
               linestyle=lineStyle[numIonTypes + i])
    ax.legend()
    ax.set_title(r"Fit {} ps, window {} $\AA$".format(fitKey, window))
    ax.set_ylabel(r"$D_{IL}(r)$  ($\AA^2$ ps$^{-1}$)")
    ax.set_xlabel(r"$r$  ($\AA$)")

