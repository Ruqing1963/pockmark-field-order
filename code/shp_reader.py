"""Minimal pure-Python ESRI shapefile (point) + dBASE (.dbf) reader.
No third-party GIS dependencies. Handles the point geometry (.shp) and
fixed-width attribute records (.dbf) of the Spatola et al. (2025) pockmark
geodatabase."""
import struct
import numpy as np


def read_points(shp_path):
    """Return (shape_type, Nx2 float array of XY) from a point shapefile."""
    b = open(shp_path, 'rb').read()
    shp_type = struct.unpack('<i', b[32:36])[0]
    pts = []; off = 100; n = len(b)
    while off < n:
        _rec, clen = struct.unpack('>ii', b[off:off+8]); off += 8
        stype = struct.unpack('<i', b[off:off+4])[0]
        if stype == 1:
            x, y = struct.unpack('<dd', b[off+4:off+20]); pts.append((x, y))
        off += clen * 2
    return shp_type, np.array(pts)


def read_dbf_column(dbf_path, field):
    """Return a list of stripped string values for one .dbf field, in record order."""
    b = open(dbf_path, 'rb').read()
    nrec = struct.unpack('<I', b[4:8])[0]
    hlen = struct.unpack('<H', b[8:10])[0]
    rlen = struct.unpack('<H', b[10:12])[0]
    pos = 32; offs = []; o = 1
    while b[pos] != 0x0D:
        nm = b[pos:pos+11].split(b'\x00')[0].decode('latin1')
        flen = b[pos+16]; offs.append((nm, flen, o)); o += flen; pos += 32
    out = []
    for i in range(nrec):
        r = b[hlen+i*rlen: hlen+(i+1)*rlen]
        for nm, flen, o2 in offs:
            if nm == field:
                out.append(r[o2:o2+flen].decode('latin1', 'replace').strip())
    return out
