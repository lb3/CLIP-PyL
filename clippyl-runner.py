#!/usr/bin/env python3.2
# -*- coding: utf-8 -*-

"""Convenience wrapper for running clippyl directly from source tree."""
#package structure emulates https://github.com/jgehrcke/python-cmdline-bootstrap
#http://gehrcke.de/2014/02/distributing-a-python-command-line-application/
#https://github.com/jgehrcke/python-cmdline-bootstrap/blob/master/bootstrap/__main__.py
#https://python-packaging-user-guide.readthedocs.org/en/latest/current.html

from clippyl.clippyl import main

if __name__ == '__main__':
    main()

