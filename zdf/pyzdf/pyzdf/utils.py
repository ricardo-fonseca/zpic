from .io import ZDFfile


def zdf_list(file_name, printRec=False):
    """
    List content of zdf file.
    """
    with ZDFfile(file_name) as zdf:
        zdf.list(printRec)


def zdf_info_grid(file_name):
    """
    Presents informations of zdf grid.
    """
    with ZDFfile(file_name) as zdf:
        if not zdf.file_type == "grid":
            raise TypeError("File is not a valid ZDF grid file")

        info = {
            "grid"      : zdf.grid_info,
            "iteration" : zdf.iteration,
        }

        return info


def zdf_read_grid(file_name):
    """
    Read grid informations
    """
    with ZDFfile(file_name) as zdf:
        if not zdf.file_type == "grid":
            raise TypeError("File is not a valid ZDF grid file")

        info = {
            "grid"      : zdf.grid_info,
            "iteration" : zdf.iteration,
        }

        data = zdf.read_dataset()

        return data, info
