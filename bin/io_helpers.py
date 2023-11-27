import gzip
import pathlib

def is_gzip(filename: pathlib.Path | str) -> bool:
    with open(filename, 'rb') as fl:
        return fl.read(2) == b'\x1f\x8b'


class FileReader:
    def __init__(self, filename: pathlib.Path | str):
        self.filename = filename
        self._file = None

    def __enter__(self):
        if is_gzip(self.filename):
            self._file = gzip.open(self.filename, 'rt')
        else:
            self._file = open(self.filename, 'rt')
        return self._file
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._file is not None:
            if not self._file.closed:
                self._file.close()


