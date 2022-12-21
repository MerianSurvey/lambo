class GaapTask():
    def __init__(self):
        from multiprocessing import Process, Manager

        self.forcedSrcCats = {}
        self.bands = list('grizy')
        manager = Manager()
        self.forcedSrcCats = manager.dict()

#     def setVal(self, cat, band, val):
#         cat[band] = val

    def run(self, band):
        self.forcedSrcCats[band] = band
#         self.forcedSrcCats[band] = band
#         return band
#         self.setVal(self.forcedSrcCats, band, band)
#         return self.forcedSrcCats

    def runAll(self, pool=None):
        """
        Run ``gaap`` photometry on all bands.

        Parameters
        ----------
        pool : multiprocessing.Pool, optional.
            The multiprocessing pool to parallelize the measurements on different bands.
        """
        if pool is not None:
            pool.map(self.run, self.bands)
            pool.close()
            pool.join()
        else:
            for band in self.bands:
                self.run(band)
