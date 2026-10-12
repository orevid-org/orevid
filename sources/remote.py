"""Reading remote parquet files via HTTP Range: only the columns that are needed.

A parquet file stores each column separately and has its index at the end of the file, so
the schema or a few columns can be read without downloading it whole.
"""
import io
import re
import urllib.request

import pyarrow.parquet as pq

HEADERS = {"User-Agent": "orevid/0.1"}


class RemoteFile(io.RawIOBase):
    def __init__(self, url):
        self.url, self.pos = url, 0
        r = urllib.request.urlopen(urllib.request.Request(url, method="HEAD", headers=HEADERS))
        self.size = int(r.headers["Content-Length"])
        self.bytes_read = 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = {0: offset, 1: self.pos + offset, 2: self.size + offset}[whence]
        return self.pos

    def read(self, n=-1):
        if n < 0:
            n = self.size - self.pos
        if n == 0 or self.pos >= self.size:
            return b""
        end = min(self.pos + n, self.size) - 1
        req = urllib.request.Request(self.url, headers={**HEADERS, "Range": f"bytes={self.pos}-{end}"})
        for attempt in range(5):
            try:
                data = urllib.request.urlopen(req, timeout=120).read()
                break
            except OSError:
                if attempt == 4:
                    raise
        self.pos += len(data)
        self.bytes_read += len(data)
        return data

    def readinto(self, b):
        d = self.read(len(b))
        b[:len(d)] = d
        return len(d)


def files(directory_url):
    html = urllib.request.urlopen(urllib.request.Request(directory_url, headers=HEADERS)).read().decode()
    return [directory_url + f for f in sorted(set(re.findall(r'href="([^"/?]+\.parquet)"', html)))]


def schema(url):
    return pq.ParquetFile(RemoteFile(url)).schema_arrow


def read(url, columns, row_filter=None):
    """Reads some columns of a remote parquet. row_filter(table) -> table, applied per row group."""
    f = RemoteFile(url)
    pf = pq.ParquetFile(f)
    parts = []
    for g in range(pf.metadata.num_row_groups):
        t = pf.read_row_group(g, columns=columns)
        parts.append(row_filter(t) if row_filter else t)
    import pyarrow as pa
    return pa.concat_tables(parts), f.bytes_read
