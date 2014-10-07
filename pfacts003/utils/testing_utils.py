"""Utilities common to Django tests"""


import httplib

def try_path(host, path, status=200):
    """Using httplib to verify a PATH is valid"""
    conn = httplib.HTTPConnection(host)
    conn.request("GET", path)

    response = conn.getresponse()

    if response.status == status:
        return True
    else:
        return False
