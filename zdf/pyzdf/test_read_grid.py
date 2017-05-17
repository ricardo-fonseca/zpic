from pyzdf.utils import zdf_read_grid

(data, info) = zdf_read_grid("J3-000500.zdf")

print(data)

print(info)
