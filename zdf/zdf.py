#!/usr/bin/python
import sys
import xdrlib
import numpy as np
import matplotlib.pyplot as plt


class ZDFfile:
    def __init__(self, file_name):
        self.__file = open(file_name, "rb")

        # Check magic number
        magic = self.__file.read(4)
        if (magic != b'ZDF0'):
            print('File is not a proper ZDF file, aborting', file=sys.stderr)
            self.__file.close

        self.__xdr = xdrlib.Unpacker("")

    def close(self):
        self.__xdr.done
        self.__file.close

    def __read_xdr_uint(self):
        self.__xdr.reset(self.__file.read(4))
        return self.__xdr.unpack_uint()

    def __read_xdr_int(self):
        self.__xdr.reset(self.__file.read(4))
        return self.__xdr.unpack_int()

    def __read_xdr_uint64(self):
        self.__xdr.reset(self.__file.read(8))
        return self.__xdr.unpack_uhyper()

    def __read_xdr_int64(self):
        self.__xdr.reset(self.__file.read(8))
        return self.__xdr.unpack_hyper()

    def __read_xdr_float(self):
        self.__xdr.reset(self.__file.read(4))
        return self.__xdr.unpack_float()

    def __read_xdr_double(self):
        self.__xdr.reset(self.__file.read(8))
        return self.__xdr.unpack_double()

    def __read_xdr_string(self):

        length = self.__read_xdr_uint()

        if (length > 0):
            # Since strings are packed as xdr_bytes there is an additional
            # length value that must be skipped
            self.__read_xdr_uint()

            # Round size up to the next multiple of 4
            size = ((length - 1) // 4 + 1) * 4
            data = self.__file.read(size)
            self.__xdr.reset(data)
            fstring = self.__xdr.unpack_fstring(length).decode()
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

    def read_record(self, skip=False):

        pos = self.__file.tell()

        # Read record id and check for EOF
        data = self.__file.read(4)
        if (data == b''):
            # If end of file return false
            return False

        # Unpack XDR record ID
        rec = dict()
        rec['pos'] = pos
        self.__xdr.reset(data)
        rec['id'] = self.__xdr.unpack_uint()
        rec['name'] = self.__read_xdr_string()
        rec['len'] = self.__read_xdr_uint64()

        if (skip):
            self.__file.seek(rec['len'], 1)

        return rec

    def read_string(self):
        rec = self.read_record()
        fstring = self.__read_xdr_string()
        return fstring

    def read_iteration(self):
        rec = self.read_record()
        iteration = dict()
        iteration['n'] = self.__read_xdr_uint()
        iteration['t'] = self.__read_xdr_double()
        iteration['tunits'] = self.__read_xdr_string()
        return iteration

    def read_grid_info(self):
        rec = self.read_record()
        info = dict()
        info['ndims'] = self.__read_xdr_uint()

        nx = []
        for i in range(info['ndims']):
            nx.append(self.__read_xdr_uint64())

        info['nx'] = nx
        info['label'] = self.__read_xdr_string()
        info['units'] = self.__read_xdr_string()
        info['has_axis'] = self.__read_xdr_uint()

        if (info['has_axis']):
            axis = []
            for i in range(info['ndims']):
                ax = dict()
                ax['type'] = self.__read_xdr_uint()
                ax['min'] = self.__read_xdr_double()
                ax['max'] = self.__read_xdr_double()
                ax['label'] = self.__read_xdr_string()
                ax['units'] = self.__read_xdr_string()
                axis.append(ax)
            info['axis'] = axis

        return info

    def read_part_info(self):
        rec = self.read_record()

        info = dict()
        info['name'] = self.__read_xdr_string()
        info['nquants'] = self.__read_xdr_uint()

        quants = []
        for i in range(info['nquants']):
            quants.append(self.__read_xdr_string())
        info['quants'] = quants

        units = dict()
        for q in quants:
            units[q] = self.__read_xdr_string()

        info['units'] = units
        info['nparts'] = self.__read_xdr_uint64()

        return info

    def read_dataset(self):
        rec = self.read_record()

        data_type = self.__read_xdr_uint()
        ndims = self.__read_xdr_uint()
        nx = []
        size = 1
        for i in range(ndims):
            nx.append(self.__read_xdr_uint64())
            size *= nx[i]

        # Read dataset and convert it to a numpy array
        # 9 - float32, 10 - float64
        if (data_type == 9):
            self.__xdr.reset(self.__file.read(size * 4))
            data = np.array(self.__xdr.unpack_farray(size,
                self.__xdr.unpack_float),
                dtype=np.float)
        elif (data_type == 10):
            self.__xdr.reset(self.__file.read(size * 8))
            data = np.array(self.__xdr.unpack_farray(size,
                self.__xdr.unpack_double),
                dtype=np.double)
        else:
            print('Unsupported datatype', file=sys.stderr)
            return False

        # Reshape dataset to the supplied dimensions
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


def zdf_info_grid(file_name):
    # Open file
    zdf = ZDFfile(file_name)
    # Check file type
    file_type = zdf.read_string()
    if (file_type != "grid"):
        print("File is not a valid ZDF grid file", file=sys.stderr)
        return False
    # Read metadata
    info = dict()
    info['grid'] = zdf.read_grid_info()
    info['iteration'] = zdf.read_iteration()
    # Close file
    zdf.close()

    return info


def zdf_read_grid(file_name):
    # Open file
    zdf = ZDFfile(file_name)
    # Check file type
    file_type = zdf.read_string()
    if (file_type != "grid"):
        print("File is not a valid ZDF grid file", file=sys.stderr)
        return False
    # Read metadata
    info = dict()
    info['grid'] = zdf.read_grid_info()
    info['iteration'] = zdf.read_iteration()
    # Read dataset
    data = zdf.read_dataset()
    # Close file
    zdf.close()

    return (data, info)


def zdf_info_particles(file_name):
    # Open file
    zdf = ZDFfile(file_name)
    # Check file type
    file_type = zdf.read_string()
    if (file_type != "particles"):
        print("File is not a valid ZDF particles file", file=sys.stderr)
        return False
    # Read metadata
    info = dict()
    info['particles'] = zdf.read_part_info()
    info['iteration'] = zdf.read_iteration()
    # Close file
    zdf.close()

    return info


def zdf_read_particles(file_name):
    # Open file
    zdf = ZDFfile(file_name)
    # Check file type
    file_type = zdf.read_string()
    if (file_type != "particles"):
        print("File is not a valid ZDF particles file", file=sys.stderr)
        return False
    # Read metadata
    info = dict()
    info['particles'] = zdf.read_part_info()
    info['iteration'] = zdf.read_iteration()
    # Read all quantities
    data = dict()
    for q in info['particles']['quants']:
        data[q] = zdf.read_dataset()
    # Close file
    zdf.close()

    return (data, info)


def zdf_list(file_name, printRec=False):
    zdf = ZDFfile(file_name)
    zdf.list(printRec)
    zdf.close()


# -----------------------------------------------------------------------------
# 2D Plotting routine
# - This is not part of this module
# -----------------------------------------------------------------------------
def plot_dataset2D(dataset, grid_info, iteration):
    fig = plt.figure()
    fig.subplots_adjust(top=0.85)
    fig.set_facecolor("#FFFFFF")

    timeLabel = r'$\sf{t = ' + str(iteration['t']) + \
        ' [' + iteration['tunits'] + r']}$'
    plotTitle = r'$\sf{' + grid_info['label'] + r'}$' + '\n' + timeLabel

    plotArea = fig.add_subplot(1, 1, 1)
    plotArea.set_title(plotTitle, fontsize=16)

    colorMap = plotArea.imshow(dataset, cmap=plt.cm.jet, 
        interpolation='nearest', origin='lower')

    colorBar = fig.colorbar(colorMap)
    colorBar.set_label(r'$\sf{' + grid_info['label'] + \
        ' [' + grid_info['units'] + r']}$', fontsize=14)

    xlabel = grid_info['axis'][0]['label'] + '[' + \
        grid_info['axis'][0]['units'] + ']'
    ylabel = grid_info['axis'][1]['label'] + '[' + \
        grid_info['axis'][1]['units'] + ']'

    plt.xlabel(r'$\sf{' + xlabel + r'}$', fontsize=14)
    plt.ylabel(r'$\sf{' + ylabel + r'}$', fontsize=14)

    return plt

# -----------------------------------------------------------------------------
# 2D Plotting routine
# - This is not part of this module
# -----------------------------------------------------------------------------


def plot_particles(dataset, info, quants):
    fig = plt.figure()
    fig.subplots_adjust(top=0.85)
    fig.set_facecolor("#FFFFFF")

    x = dataset[quants[1]]
    y = dataset[quants[0]]

    plt.plot(x, y, 'r.', alpha=0.5)

    t = str(info["iteration"]["t"])
    tunits = str(info["iteration"]["tunits"])

    title = info['particles']['name'] + '  ' + quants[1] + ' ' + quants[0]

    timeLabel = r'$\sf{t = ' + t + ' [' + tunits + r']}$'
    plt.title(r'$\sf{' + title + r'}$' + '\n' + timeLabel)

    xlabel = quants[1] + '[' + \
        info['particles']['units'][quants[1]] + ']'
    ylabel = quants[0] + '[' + \
        info['particles']['units'][quants[0]] + ']'

    plt.xlabel(r'$\sf{' + xlabel + r'}$', fontsize=14)
    plt.ylabel(r'$\sf{' + ylabel + r'}$', fontsize=14)

    return plt


# -----------------------------------------------------------------------------
# Main code
# -----------------------------------------------------------------------------

# List file contents
# zdf_list( "J3-000500.zdf", printRec = True )

# Print grid metadata
# print( zdf_info_grid( "J3-000500.zdf" ) )

# Plot grid


# (data, info) = zdf_read_grid("J3-000500.zdf")
# plt = plot_dataset2D(data, info['grid'], info['iteration'])
# plt.show()

# zdf_list("particles-electrons-001200.zdf", printRec = True)
# (particles, info) = zdf_read_particles("particles-electrons-001200.zdf")
# plt2 = plot_particles(particles, info, ('u1', 'x1'))
# plt2.show()
