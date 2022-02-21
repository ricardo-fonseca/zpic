"""
Python module for reading ZDF data files

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

import sys
import numpy as np

class ZDF_Record:
    """ZDF_Record()

    ZDF Datafile record information

    Attributes
    ----------
    pos : int64
        Position in ZDF file
    id : int
        Type id
    name : str
        Record name
    len : uint64
        Additional record length in bytes
    """
    def __init__( self, pos, id, name, len ):
        self.pos   = pos
        self.id    = id
        self.name  = name
        self.len   = len
    
    def version( self ):
        """version()

        Gets record version number

        Returns
        -------
        version : int
            Record version number
        """
        return self.id & 0x0000FFFF

    def type( self ):
        """type( typeTag )

        Gets record type from tag

        Returns
        -------
        type : {'int', 'double', 'string', 'dataset', 'iteration', 'grid_info', 'part_info','unknown'}
            Type of record
        """

        typeID = self.id & 0xFFFF0000

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


class ZDF_Iteration:
    """ZDF_Iteration()
    
    Class describing iteration information.

    Attributes
    ----------
    n : int
        Iteration value
    t : float
        Time value
    tunits : str
        Units used for time value
    """
    def __init__( self ):
        self.n = 0
        self.t = 0.0
        self.tunits = ""

class ZDF_Grid_Axis:
    """ZDF_Grid_Axis()

    Class describing grid axis

    Attributes
    ----------
    type : {0,1,2}
        Axis type, must be one of 0 (linear), 1 (log10) or 2 (log2)
    min : float
        Minimum axis value
    max : float
        Maximum axis value
    label : str
        Axis label
    units : str
        Axis values data units
    """
    def __init__( self ):
        self.type = 0
        self.min = 0.
        self.max = 0.
        self.label = ''
        self.units = ''

class ZDF_Grid_Info:
    """ZDF_Grid_Info()

    Grid dataset information

    Attributes
    ----------
    ndims : int
        Dimensionality of dataset
    nx : list of int (ndims)
        Number of grid points in each direction
    label : str
        Dataset label
    units : str
        Dataset units
    has_axis : bool
        True if the dataset includes axis information
    axis : list of ZDF_Grid_Axis (ndims)
        Information on each axis, when available
    """
    def __init__( self ):
        self.ndims = 0
        self.nx = []
        self.label = ''
        self.units = ''
        self.has_axis = 0
        self.axis = []

class ZDF_Part_Info:
    """ZDF_Part_Info()

    Particle dataset information

    Attributes
    ----------
    name : str
        Particle dataset name
    nquants : int
        Number of quantities per particle
    quants : list of str (nquants)
        Name of individual quantities
    units : dictionary
        Units for each quantity
    nparts: int
        Number of particles in dataset
    """
    def __init__( self ):
        self.name = ''
        self.nquants = 0
        self.quants = []
        self.units = dict()
        self.nparts = 0

class ZDFfile:
    """ZDFfile( file_name )
    
    ZDF data file class

    Parameters
    ----------
    file_name : str
        File name of ZDF data file, should include path
    """
    def __init__(self, file_name):
        self.__file = open(file_name, "rb")

        # Check magic number
        magic = self.__file.read(4)
        if (magic != b'ZDF1'):
            print('File is not a proper ZDF file, aborting', file=sys.stderr)
            self.__file.close

    def close(self):
        """close()

        Closes ZDF file
        """
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

# -----------------------------------------------------------------------------
# Low level interfaces
# -----------------------------------------------------------------------------

    def read_record(self, skip=False):
        """read_record(skip=False)

        Reads current record information from file

        Parameters
        ----------
        skip : bool, optional
            If set to True, skip to next record after reading record
            header data
        """

        pos = self.__file.tell()

        # Read record id and check for EOF
        id = self.__read_uint32()

        if (id is False):
            # If end of file return false
            return False

        # Read record ID
        name = self.__read_string()
        len = self.__read_uint64()
        rec = ZDF_Record( pos, id, name, len )

        # If requested, skip over to next record
        if (skip):
            self.__file.seek(rec.len, 1)

        return rec

    def read_string(self):
        """read_string()

        Reads string record from data file

        Returns
        -------
        string : str
            String data
        """
        rec = self.read_record()
        fstring = self.__read_string()
        return fstring

    def read_iteration(self):
        """read_iteration()

        Read iteration record from data file

        Returns
        -------
        iteration : ZDF_Iteration()
            Iteration data
        """
        rec = self.read_record()
        iteration = ZDF_Iteration()
        iteration.n = self.__read_int32()
        iteration.t = self.__read_float64()
        iteration.tunits = self.__read_string()
        return iteration

    def read_grid_info(self):
        """read_grid_info()

        Read grid information record from data file

        Returns
        -------
        iteration : ZDF_Grid_Info()
            Grid information data
        """
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
        """read_part_info()

        Read particle information record from data file

        Returns
        -------
        iteration : ZDF_Part_Info()
            Iteration data
        """
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
        """read_dataset()

        Read dataset from data file

        Returns
        -------
        data : numpy.ndarray
            Numpy ndarray with data
        """
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
        """list( printRec=False )

        Gets a list of file contents and optionally prints it to screen

        Parameters
        ----------
        printRec : bool, optional
            If set to True will print all records found in the file,
            defaults to False.
        """
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
                    rec.pos, rec.len,
                    self.record_type(rec.id), rec.name))

        return rec_list

# -----------------------------------------------------------------------------
# High level interfaces
# -----------------------------------------------------------------------------

class ZDF_Info:
    """ZDF_Info()
    ZDF File information

    Attributes
    ----------
    type : {'grid','particles'}
        Type of ZDF file
    grid : ZDF_Grid_Info
        Grid information for grid files
    particles : ZDF_Part_Info
        Particle information for particle files
    iteration : ZDF_Iteration
        Iteration information

    """
    def __init__(self):
        self.type = ""
        self.grid = None
        self.particles = None
        self.iteration = None


def info( file_name ):
    """info( file_name )

    Gets metadata for a ZDF file

    Parameters
    ----------
    file_name : str
        File name of ZDF data file, should include path
    
    Returns
    -------
    info : ZDF_Info
        File information. If file is invalid False is returned.
    """
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
    """read( file_name )

    Reads all data in a ZDF file

    Parameters
    ----------
    file_name : str
        File name of ZDF data file, should include path
    
    Returns
    -------
    (data, info) : ( numpy.ndarray | dictionary, ZDF_Info )
        Tuple containing file data and metadata. Data will be a
        numpy.ndarray for grid data, and a dictionary of numpy.array for
        particle data (one entry per quantity). Metadata is returned as a
        ZDF_Info object. If file is invalid False is returned.
    """
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
    """list( printRec=False )

    Gets a list of file contents and optionally print it to screen

    Parameters
    ----------
    file_name : str
        File name of ZDF data file, should include path
        
    printRec : bool, optional
        If set to True will print all records found in the file,
        defaults to False.
    """
    zdf = ZDFfile(file_name)
    zdf.list(printRec)
    zdf.close()

