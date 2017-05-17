import xdrlib
import numpy as np

# defines global imports (needed for * imports)
__all__ = ["ZDFfile"]

# constants - byte lengths for zdf content
MAGIC_NUMBER = 4
LEN_32 = 4
LEN_64 = 8


class _XDRstream:
    _xdr = xdrlib.Unpacker("")

    def _read_xdr_uint(self):
        self._xdr.reset(self._file.read(LEN_32))
        return self._xdr.unpack_uint()

    def _read_xdr_int(self):
        self._xdr.reset(self._file.read(LEN_32))
        return self._xdr.unpack_int()

    def _read_xdr_uint64(self):
        self._xdr.reset(self._file.read(LEN_64))
        return self._xdr.unpack_uhyper()

    def _read_xdr_int64(self):
        self._xdr.reset(self._file.read(LEN_64))
        return self._xdr.unpack_hyper()

    def _read_xdr_float(self):
        self._xdr.reset(self._file.read(LEN_32))
        return self._xdr.unpack_float()

    def _read_xdr_double(self):
        self._xdr.reset(self._file.read(LEN_64))
        return self._xdr.unpack_double()

    def _read_xdr_string(self):
        length = self._read_xdr_uint()

        if (length > 0):
            # Since strings are packed as xdr_bytes there is an additional
            #  length value that must be skipped
            self._read_xdr_uint()

            # Round size up to the next multiple of LEN_32
            size = ((length - 1) // LEN_32 + 1) * LEN_32
            data = self._file.read(size)
            self._xdr.reset(data)
            fstring = self._xdr.unpack_fstring(length).decode()
        else:
            fstring = ""

        return fstring


class _FileContextManager:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is not None:
            self.close
            return True

        return False


class ZDFfile(_XDRstream, _FileContextManager):
    _types = {
        0x00010000 : "int",
        0x00020000 : "double",
        0x00030000 : "string",
        0x00100000 : "dataset",
        0x00200000 : "iteration",
        0x00210000 : "grid_info",
    }
    _file_type = None

    def __init__(self, file_name):
        self._file = open(file_name, "rb")

        # check magic number
        magic = self._file.read(MAGIC_NUMBER)

        if str(magic, "utf-8") != "ZDF0":
            self._file.close
            raise IOError("File is not a proper ZDF file ... aborting")

    @property
    def close(self):
        self._xdr.done
        self._file.close
        return True

    def _record_type(self, typeTag):
        version = typeTag & 0x0000FFFF
        typeID  = typeTag & 0xFFFF0000

        if typeID in self._types:
            return self._types[typeID]
        else:
            return "unknown"

    @property
    def file_type(self):
        if not self._file_type:
            self._file_type = self.read_string()
            return self._file_type
        else:
            return self._file_type

    @property
    def grid_info(self):
        rec = self.read_record()
        info = dict()
        info['ndims'] = self._read_xdr_uint()

        nx = []
        for i in range(info['ndims']):
            nx.append(self._read_xdr_uint64())

        info['nx'] = nx
        info['label'] = self._read_xdr_string()
        info['units'] = self._read_xdr_string()
        info['has_axis'] = self._read_xdr_uint()

        if info['has_axis']:
            axis = []
            for i in range(info['ndims']):
                ax = dict()
                ax['type'] = self._read_xdr_uint()
                ax['min'] = self._read_xdr_double()
                ax['max'] = self._read_xdr_double()
                ax['label'] = self._read_xdr_string()
                ax['units'] = self._read_xdr_string()
                axis.append(ax)
            info['axis'] = axis

        return info

    @property
    def iteration(self):
        rec = self.read_record()
        iteration = dict()
        iteration['n'] = self._read_xdr_uint()
        iteration['t'] = self._read_xdr_double()
        iteration['tunits'] = self._read_xdr_string()
        return iteration

    def read_string(self):
        rec     = self.read_record()
        fstring = self._read_xdr_string()
        return fstring

    def read_record(self, skip=False):
        pos = self._file.tell()

        # Read record id and check for EOF
        data = self._file.read(4)
        if data == b'':
            # If end of file return false
            return False

        # Unpack XDR record ID
        rec = dict()
        rec['pos']  = pos
        self._xdr.reset(data)

        rec['id']   = self._xdr.unpack_uint()
        rec['name'] = self._read_xdr_string()
        rec['len']  = self._read_xdr_uint64()

        if skip:
            self._file.seek(rec['len'], 1)

        return rec

    def read_records(self, skip=False):
        while True:
            rec = self.read_record(skip)
            if not rec:
                raise StopIteration
            else:
                yield rec

    def list(self, printRec=False):
        # Position file after magic number
        self._file.seek(MAGIC_NUMBER)

        # header for list option
        print(
            "{:^16s}".format("Position") + "|" +
            "{:^16s}".format("Size (bytes)") + "|" +
            "{:^16s}".format("Type") + "|" +
            "{:^16s}".format("Name") + "\n" +
            "-" * (16 * 4 + 3)
        )

        # iterate over records - skip content
        for rec in self.read_records(skip=True):
            rec["id"] = self._record_type(rec["id"])
            print(" {pos:#014x} | {len:#014x} |".format(**rec) +
                  " {id:^14s} | {name:^14s}".format(**rec))

    def read_dataset(self):
        rec = self.read_record()

        data_type = self._read_xdr_uint()
        ndims     = self._read_xdr_uint()
        nx = []
        size = 1
        for i in range(ndims):
            nx.append(self._read_xdr_uint64())
            size *= nx[i]

        # Read dataset and convert it to a numpy array
        # 9 - float32, 10 - float64
        if data_type == 9:
            self._xdr.reset(self._file.read(size * LEN_32))
            data = np.array(
                self._xdr.unpack_farray(size, self._xdr.unpack_float),
                dtype=np.float
            )

        elif data_type == 10:
            self._xdr.reset(self._file.read(size * LEN_64))
            data = np.array(
                self._xdr.unpack_farray(size, self._xdr.unpack_double),
                dtype=np.double
            )

        else:
            raise TypeError("Unsupported datatype")

        # Reshape dataset to the supplied dimensions
        data.shape = nx

        return data
