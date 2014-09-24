sdDInfty2 = {}
aveW = np.arange(0.01, 0.51, dr) #nm
for fit in sdD:
  sdDInfty2[fit] = []
  for i in range(sdD[fit].shape[0]):
    sdDInfty2[fit].append([np.mean(sdD[fit][i, halfCellIndex - aW / dr: halfCellIndex + aW / dr]) for aW in aveW])

for fit in sdDInfty2:
  plt.figure()
  plt.title(fit)
  if (args.color is not None):
    plt.gca().set_color_cycle(args.color[numIonTypes:])
  for i, sdDInf2 in enumerate(sdDInfty2[fit]):
    plt.plot(aveW * Const.nm2AA, sdDInf2, linestyle=lineStyle[numIonTypes + i], label=label[numIonTypes + i])
  plt.legend()
  plt.xlabel(r"average window  ($\AA$)")
  plt.ylabel(r"$D_{IL}(\infty)$  ($\AA^2$ $ps^{-1}$)")
