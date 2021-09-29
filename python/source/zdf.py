"""
Copyright (C) 2017 Instituto Superior Tecnico

This file is part of the ZPIC Educational code suite

The ZPIC Educational code suite is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The ZPIC Educational code suite is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with the ZPIC Educational code suite. If not, see <http://www.gnu.org/licenses/>.
"""

#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

class ZDF_Iteration:
    def __init__( self ):
        self.n = 0
        self.t = 0.
        self.tunits = ""

class ZDF_Grid_Axis:
    def __init__( self ):
        self.type = 0
        self.min = 0.
        self.max = 0.
        self.label = ''
        self.units = ''

class ZDF_Grid_Info:
    def __init__( self ):
        self.ndims = 0
        self.nx = []
        self.label = ''
        self.units = ''
        self.has_axis = 0
        self.axis = []

class ZDF_Part_Info:
    def __init__( self ):
        self.name = ''
        self.nquants = 0
        self.quants = []
        self.units = dict()
        self.nparts = 0

class ZDFfile:
    def __init__(self, file_name):
        self.__file = open(file_name, "rb")

        # Check magic number
        magic = self.__file.read(4)
        if (magic != b'ZDF1'):
            print('File is not a proper ZDF file, aborting', file=sys.stderr)
            self.__file.close

    def close(self):
        self.__file.close

    def __read_uint32(self):
        data = np.fromfile(self.__file,dtype='<u4',count=1)
        if ( data.size == 0 ):
            return False
        else:
            return data[0]

    def __read_int32(self):
        return np.fromfile(self.__file,dtype='<i4',count=1)[0]

    def __read_uint64(self):
        return np.fromfile(self.__file,dtype='<u8',count=1)[0]

    def __read_int64(self):
        return np.fromfile(self.__file,dtype='<i8',count=1)[0]

    def __read_float32(self):
        return np.fromfile(self.__file,dtype='<f4',count=1)[0]

    def __read_float64(self):
        return np.fromfile(self.__file,dtype='<f8',count=1)[0]

    def __read_string(self):

        length = self.__read_uint32()

        if (length > 0):
            data = self.__file.read(length)
            fstring = data.decode()

            # Data is always written in blocks of 4 byt
            pad = ((length - 1) // 4 + 1) * 4 - length
            if ( pad > 0 ):
                self.__file.seek(pad,1)
        else:
            fstring = ""

        return fstring

    def record_type(self, typeTag):
        version = typeTag & 0x0000FFFF
        typeID = typeTag & 0xFFFF0000

        types = {0x00010000: "int",
                 0x00020000: "double",
                 0x00030000: "string",
                 0x00100000: "dataset",
                 0x00200000: "iteration",
                 0x00210000: "grid_info",
                 0x00220000: "part_info", }

        if (typeID in types):
            return types[typeID]
        else:
            return "unknown"

# -----------------------------------------------------------------------------
# Low level interfaces
# -----------------------------------------------------------------------------

    def read_record(self, skip=False):

        pos = self.__file.tell()

        # Read record id and check for EOF
        id = self.__read_uint32()

        if (id is False):
            # If end of file return false
            return False

        # Read record ID
        rec = dict()
        rec['pos'] = pos
        rec['id']  = id
        rec['name'] = self.__read_string()
        rec['len']  = self.__read_uint64()

        if (skip):
            self.__file.seek(rec['len'], 1)

        return rec

    def read_string(self):
        rec = self.read_record()
        fstring = self.__read_string()
        return fstring

    def read_iteration(self):
        rec = self.read_record()
        iteration = ZDF_Iteration()
        iteration.n = self.__read_int32()
        iteration.t = self.__read_float64()
        iteration.tunits = self.__read_string()
        return iteration

    def read_grid_info(self):
        rec = self.read_record()
        info = ZDF_Grid_Info()
        info.ndims = self.__read_uint32()

        for i in range(info.ndims):
            info.nx.append(self.__read_uint64())

        info.label = self.__read_string()
        info.units = self.__read_string()
        info.has_axis = self.__read_int32()

        if ( info.has_axis ):
            for i in range(info.ndims):
                ax = ZDF_Grid_Axis()
                ax.type  = self.__read_int32()
                ax.min   = self.__read_float64()
                ax.max   = self.__read_float64()
                ax.label = self.__read_string()
                ax.units = self.__read_string()
                info.axis.append(ax)

        return info

    def read_part_info(self):
        rec = self.read_record()

        info = ZDF_Part_Info()
        info.name    = self.__read_string()
        info.nquants = self.__read_uint32()

        for i in range(info.nquants):
            info.quants.append(self.__read_string())

        for q in info.quants:
            info.units[q] = self.__read_string()

        info.nparts = self.__read_uint64()

        return info

    def read_dataset(self):
        rec = self.read_record()

        data_type = self.__read_int32()
        ndims = self.__read_uint32()
        nx = []
        size = np.uint64(1)

        for i in range(ndims):
            nx.append(self.__read_uint64())
            size *= nx[i]

        # Read dataset and convert it to a numpy array
        # 9 - float32, 10 - float64
        if (data_type == 9):
            data = np.fromfile(self.__file,dtype='float32',count=size)
        elif (data_type == 10):
            data = np.fromfile(self.__file,dtype='float64',count=size)
        else:
            print('Unsupported datatype', file=sys.stderr)
            return False

        # Reshape dataset to the supplied dimensions
        nx.reverse()
        data.shape = nx

        return data

    def list(self, printRec=False):
        # Position file after magic number
        self.__file.seek(4)

        rec_list = []
        while True:
            rec = self.read_record(skip=True)
            if (rec is False):
                break
            else:
                rec_list.append(rec)

        if (printRec and (len(rec_list) > 0)):
            print('Position     Size(bytes)  Type        Name')
            print('-------------------------------------------------')
            for rec in rec_list:
                print('{:#010x}   {:#010x}   {:10}  {}'.format(
                    rec['pos'], rec['len'],
                    self.record_type(rec['id']), rec['name']))

        return rec_list

# -----------------------------------------------------------------------------
# High level interfaces
# -----------------------------------------------------------------------------

class ZDF_Info:
    """ZDF File information"""
    def __init__(self):
        self.type = ""
        self.grid = None
        self.particles = None
        self.iteration = None


def info( file_name ):
    # Open file
    zdf = ZDFfile( file_name )

    # Check file type and read metadata
    info = ZDF_Info()
    info.type = zdf.read_string()
    if ( info.type == "grid" ):
        info.grid = zdf.read_grid_info()
    elif ( info.type == "particles" ):
        info.particles = zdf.read_part_info()
    else:
        print("File is not a valid ZDF grid or particles file", file=sys.stderr)
        zdf.close()
        return False

    # Read iteration info
    info.iteration = zdf.read_iteration()

    # Close file
    zdf.close()

    return info

def read( file_name ):
    # Open file
    zdf = ZDFfile( file_name )

    # Check file type and read metadata
    info = ZDF_Info()

    info.type = zdf.read_string()
    if ( info.type == "grid" ):
        info.grid = zdf.read_grid_info()
        info.iteration = zdf.read_iteration()
        data = zdf.read_dataset()
    elif ( info.type == "particles" ):
        info.particles = zdf.read_part_info()
        info.iteration = zdf.read_iteration()

        # Read all quantities
        data = dict()
        for q in info.particles.quants:
            data[q] = zdf.read_dataset()
    else:
        print("File is not a valid ZDF grid or particles file", file=sys.stderr)
        zdf.close()
        return False

    #Close file
    zdf.close()

    return (data,info)

def list(file_name, printRec=False):
    zdf = ZDFfile(file_name)
    zdf.list(printRec)
    zdf.close()

