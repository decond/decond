import h5py
import os.path
from ._version import __version__


class DecondFile(h5py.File):
    """
    A container file of all the analyzed data
    """
    def __init__(self, name, mode='a', **kwarg):
        isFileExistAtFirst = os.path.isfile(name)
        super().__init__(name, mode, **kwarg)
        if isFileExistAtFirst and mode in ('r', 'r+', 'a'):
            self._checkVersion()
        if ((not isFileExistAtFirst and mode is 'a') or
                mode in ('w', 'w-', 'x')):
            self._initEmptyData()

    def _checkVersion(self):
        self.version = DecondFile._readAttrStr(self, 'version')
        if self.version is not None:
            (fmajor, fminor, fpatch) = (self.version.decode().split(sep='.'))
            (MAJOR, MINOR, PATCH) = (__version__.split(sep='.'))
            if fmajor < MAJOR:
                raise DecondFile.IOError(
                        "File " + self.filename +
                        " is of version " + fmajor + ".X.X, " +
                        "while this program requires at least " +
                        MAJOR + ".X.X")
        else:
            raise DecondFile.IOError(
                    "File " + self.filename +
                    " has no version number. " +
                    "This program requires files of at least " +
                    "version " + MAJOR + ".X.X")

    def _initEmptyData(self):
        self['numSample'] = 0

    @staticmethod
    def _readDset(base, dset):
        return base[dset][...] if dset in base else None

    @staticmethod
    def _readAttrStr(base, attr):
        return base.attrs[attr] if attr in base.attrs else None

    @staticmethod
    def _readGroupToDict(base, group):
        D = {}
        if group in base:
            for k, v in group.items():
                D[k] = v[...]
            return D
        else:
            return None

    def addSample(self, samples):
        if type(samples) is not list:
            samples = [samples]

        if self['numSample'][...] == 0:
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

    def writeBuffer(self):
        self.attrs['version'] = self.buffer.version
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
        self['nD'] = self.buffer.nD
        self['nD_err'] = self.buffer.nD_err
        self['nDCesaroTotal'] = self.buffer.nDCesaroTotal
        self['nDCesaroTotal_err'] = self.buffer.nDCesaroTotal_err
        self['nDTotal'] = self.buffer.nDTotal
        self['nDTotal_err'] = self.buffer.nDTotal_err

        if self.buffer.spatialDec is not None:
            self.buffer.spatialDec.write(self, 'spatialDec')

        if self.buffer.energyDec is not None:
            self.buffer.energyDec.write(self, 'energyDec')

    class Buffer:
        def __init__(self, file):
            if type(file) is DecondFile:
                self._file = file
                closeFilePlease = False
            else:
                self._file = h5py.File(file, 'r')
                closeFilePlease = True

            self.version = self._file.attrs['version']
            self.numSample = DecondFile._readDset(self._file, 'numSample')
            self.charge = self._file['charge'][...]
            self.numMol = self._file['numMol'][...]
            self.volume = self._file['volume'][...]
            self.volume_unit = self._file['volume'].attrs['unit']
            self.volume_err = DecondFile._readDset(self._file, 'volume_err')
            self.timeLags = self._file['timeLags'][...]

            self.nCorr = self._file['nCorr'][...]
            self.nCorr_err = DecondFile._readDset(self._file, 'nCorr_err')

            self.nDCesaro = DecondFile._readDset(self._file, 'nDCesaro')
            self.nDCesaro_err = DecondFile._readDset(
                    self._file, 'nDCesaro_err')

            self.nD = DecondFile._readGroupToDict(self._file, 'nD')
            self.nD_err = DecondFile._readGroupToDict(self._file, 'nD_err')

            self.nDCesaroTotal = DecondFile._readGroupToDict(
                    self._file, 'nDCesaroTotal')
            self.nDCesaroTotal_err = DecondFile._readGroupToDict(
                    self._file, 'nDCesaroTotal_err')

            self.nDTotal = DecondFile._readGroupToDict(self._file, 'nDTotal')
            self.nDTotal_err = DecondFile._readGroupToDict(
                    self._file, 'nDTotal_err')

            self.spatialDec = DecondFile.DecData(self._file, 'spatialDec')

            self.energyDec = DecondFile.DecData(self._file, 'energyDec')

            if closeFilePlease:
                self._file.close()

    class DecData:
        def __init__(self, base, group):
            if group in base:
                decGroup = base['group']
            else:
                return None

            self.decBins = decGroup['decBins'][:]
            self.decCorr = decGroup['decCorr'][:, :, :]
            self.decPairCount = decGroup['decPairCount'][:, :]

            self.decBins_unit = DecondFile._readAttrStr(
                    decGroup['decBins'], 'unit')

            self.decCorr_unit = DecondFile._readAttrStr(
                    decGroup['decCorr'], 'unit')

            self.decCorr_err = DecondFile._readDset(decGroup, 'decCorr_err')

            self.decPairCount_err = DecondFile._readDset(
                    decGroup, 'decPairCount_err')

            self.decDCesaro = DecondFile._readGroupToDict(
                    decGroup, 'decDCesaro')
            self.decDCesaro_err = DecondFile._readGroupToDict(
                    decGroup, 'decDCesaro_rr')
            if self.decDCesaro is not None:
                self.decDCesaro_unit = DecondFile._readAttrStr(
                        decGroup['decDCesaro'], 'unit')

            self.decD = DecondFile._readGroupToDict(decGroup, 'decD')
            self.decD_err = DecondFile._readGroupToDict(decGroup, 'decD_err')

            if self.decD is not None:
                self.decD_unit = DecondFile._readAttrStr(
                        decGroup['decD'], 'unit')

        def write(self, base, group):
            decGroup = base.require_group(group)
            decGroup['decBins'] = self.decBins
            decGroup['decCorr'] = self.decCorr
            decGroup['decPairCount'] = self.decPairCount

            decGroup['decBins'].attrs['unit'] = self.decBins_unit
            decGroup['decCorr'].attrs['unit'] = self.decCorr_unit
            decGroup['decCorr_err'] = self.decCorr_err

            decGroup['decPairCount_err'] = self.decPairCount_err

            decGroup['decDCesaro'] = self.decDCesaro
            decGroup['decDCesaro_err'] = self.decDCesaro_err
            decGroup['decDCesaro'].attrs['unit'] = self.decDCesaro_unit

            decGroup['decD'] = self.decD
            decGroup['decD_err'] = self.decD_err
            decGroup['decD'].attrs['unit'] = self.decD_unit

    class IOError(Exception):
        pass
