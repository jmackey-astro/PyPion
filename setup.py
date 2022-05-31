#!/usr/bin/env python

from setuptools import setup, __version__
from setuptools.command.install import install
from pkg_resources import parse_version
import subprocess
import sys
import os

class CustomInstall(install):
    def run(self):
        command = "mkdir ~/.silo_test"
        process = subprocess.Popen(command, shell=True)
        process.wait()
        install.run(self)

minimum_version = parse_version('30.4.0')

if parse_version(__version__) < minimum_version:
    raise RuntimeError("Package setuptools must be at least version {}".format(minimum_version))

setup(cmdclass={'install': CustomInstall})

#os.system('mkdir ~/.silo_test')
