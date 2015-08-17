import h5py
import os.path
import numpy as np
from ._version import __version__


class DecondIOError(Exception):
    pass


class DecData:
    def __init__(self, decGroup):
        self.decBins = decGroup['decBins'][:]
        self.decCorr = decGroup['decCorr'][:, :, :]
        self.decPairCount = decGroup['decPairCount'][:, :]

        self.decBins_unit = DecondFile.readAttrStr(decGroup['decBins'],
                                                   'unit')

        self.decCorr_unit = DecondFile.readAttrStr(decGroup['decCorr'],
                                                   'unit')

        self.decCorr_err = DecondFile.readOptDset(decGroup, 'decCorr_err')

        self.decPairCount_err = DecondFile.readOptDset(decGroup,
                                                       'decPairCount_err')

        self.decDCesaro = DecondFile.readGroupToDict(decGroup,
                                                     'decDCesaro')
        self.decDCesaro_err = DecondFile.readGroupToDict(decGroup,
                                                         'decDCesaro_rr')
        self.decDCesaro_unit = DecondFile.readAttrStr(decGroup['decDCesaro'],
                                                      'unit')

        self.decD = DecondFile.readGroupToDict(decGroup, 'decD')
        self.decD_err = DecondFile.readGroupToDict(decGroup, 'decD_err')
        self.decD_unit = DecondFile.readAttrStr(decGroup['decD'], 'unit')


class DecondFile(h5py.File):
    """
    A container file of all the analyzed data
    """
    def __init__(self, name, mode, **kwarg):
        isFileExistAtFirst = os.path.isfile(name)
        super().__init__(name, mode, **kwarg)
        if isFileExistAtFirst and mode in ('r', 'r+', 'a'):
            self._checkVersion()
        if ((not isFileExistAtFirst and mode is 'a') or
                mode in ('w', 'w-', 'x')):
            self._initEmptyData()

    def _checkVersion(self):
        self.version = DecondFile.readAttrStr(self, 'version')
        if self.version is not None:
            (fmajor, fminor, fpatch) = (self.version.decode().split(sep='.'))
            (MAJOR, MINOR, PATCH) = (__version__.split(sep='.'))
            if fmajor < MAJOR:
                raise DecondIOError("File " + self.filename +
                                    " is of version " +
                                    fmajor + ".X.X, " +
                                    "while this program requires at least " +
                                    MAJOR + ".X.X")
        else:
            raise DecondIOError("File " + self.filename +
                                " has no version number. " +
                                "This program requires files of at least " +
                                "version " + MAJOR + ".X.X")

    def _initEmptyData(self):
        self.attrs['version'] = np.string_(__version__)
        self['numSample'] = 0

    @staticmethod
    def readOptDset(base, dset):
        return base[dset][...] if dset in base else None

    @staticmethod
    def readAttrStr(base, attr):
        return base.attrs[attr] if attr in base.attrs else None

    @staticmethod
    def readGroupToDict(base, group):
        D = {}
        if group in base:
            for k, v in group.items():
                D[k] = v[...]
            return D
        else:
            return None

    def addSample(self, samples):
        if self['numSample'] == 0:
            # new file, read in the first sample
            with DecondFile(samples[0]) as f:
                f.readBuffer()
                self.buffer = f.buffer
            begin = 1
        else:
            # read data from the existing decond file
            self.readBuffer()
            begin = 0

        for sample in samples[begin:]:
            with h5py.File(sample, 'r') as f:
                pass

        self.writeBuffer()

    def readBuffer(self):
        self.buffer = DecondFile.Buffer(self)

    @todo
    def writeBuffer(self):
        pass

    class Buffer:
        def __init__(self, file):
            if type(file) is DecondFile:
                self._file = file
                closeFilePlease = False
            else:
                self._file = h5py.File(file, 'r')
                closeFilePlease = True

            self.version = self._file.attrs['version']
            self.numSample = DecondFile.readOptDset(self._file, 'numSample')
            self.charge = self._file['charge']
            self.numMol = self._file['numMol']
            self.volume = self._file['volume']
            self.volume_unit = self._file['volume'].attrs['unit']
            self.volume_err = DecondFile.readOptDset(self._file, 'volume_err')
            self.timeLags = self._file['timeLags']

            self.nCorr = self.file['nCorr']
            self.nCorr_err = DecondFile.readOptDset(self._file, 'nCorr_err')

            self.nDCesaro = DecondFile.readOptDset(self._file, 'nDCesaro')
            self.nDCesaro_err = DecondFile.readOptDset(self._file,
                                                       'nDCesaro_err')
            self.nD = DecondFile.readGroupToDict(self._file, 'nD')
            self.nD_err = DecondFile.readGroupToDict(self._file, 'nD_err')

            self.nDCesaroTotal = DecondFile.readGroupToDict(self._file,
                                                            'nDCesaroTotal')
            self.nDCesaroTotal_err = DecondFile.readGroupToDict(self._file,
                                                        'nDCesaroTotal_err')

            self.nDTotal = DecondFile.readGroupToDict(self._file, 'nDTotal')
            self.nDTotal_err = DecondFile.readGroupToDict(self._file,
                                                          'nDTotal_err')

            self.spatialDec = DecData(self._file['spatialDec'])
            self.energyDec = DecData(self._file['energyDec'])

            if closeFilePlease:
                self._file.close()
