''' This uses the imexam package and makes suitable for use in this package'''

from imexam import connect as Connect, ds9_viewer
from imexam.imexamine import Imexamine
from imexam.util import set_logging
import os
import shutil
import sys
import numpy as np
from subprocess import Popen
import time
import logging
from tempfile import mkdtemp
import warnings

try:
    from imexam import imexamxpa
    have_xpa = True
    from imexam.ds9_viewer import ds9 as DS9
except ImportError:
    have_xpa = False

class connect(Connect):

    def __init__(self, target=None, path=None, viewer="ds9",
                 wait_time=10, quit_window=True, port=None):
        """Initialize the imexam control object."""
        _possible_viewers = []

        if have_xpa:
            _possible_viewers.append("ds9")

        self._viewer = viewer.lower()

        if (self._viewer not in _possible_viewers or
                len(_possible_viewers) == 0):
            warnings.warn("**Unsupported viewer, check your installed "
                          "packages**\n")
            raise NotImplementedError

        # init sets empty data array until we can load or check viewer
        self.exam = Imexamine()

        if 'ds9' in self._viewer:
            self.window = ds9(
                target=target,
                path=path,
                wait_time=wait_time,
                quit_ds9_on_del=quit_window)
            self._event_driven_exam = False  # use the imexam loop


            # alter the exam.imexam_option_funcs{} here through the viewer code
            # if you want to change key+function associations
            # self.window._reassign_keys(imexam_dict)

        self.logfile = 'imexam_log.txt'  # default logfile name
        self.log = set_logging()  # points to the package logger
        self._current_slice = None
        self._current_frame = None

class ds9(DS9):

    def _run_unixonly_ds9(self):
        """start new ds9 window and connect to object using a unix socket.

        Notes
        -----
        When the xpa method in libxpa parses a given template as a unix
        socket, it checks if the template string starts with tmpdir
        (from env["XPA_TMPDIR"] or default to /tmp/.xpa). This can make
        having multiple instances of ds9 a bit difficult, but if you give it
        unique names or use the inet address you should be fine

        For unix only, we run ds9 with XPA_TMPDIR set to temporary directory
        whose prefix start with /tmp/xpa (eg, /tmp/xpa_sf23f), them set
        os.environ["XPA_TMPDIR"] (which affects xpa set and/or get command
        from python) to /tmp/xpa.
        """
        env = os.environ
        wait_time = self.wait_time

        self._tmpd_name = mkdtemp(
            prefix="xpa_" +
            env.get(
                "USER",
                ""),
            dir="/tmp")

        # this is the first directory the servers looks for on the path
        env["XPA_TMPDIR"] = self._tmpd_name

        unix_name = "{0:s}/.IMT".format(self._tmpd_name)

        # that should be unique enough, something better?
        title = str(time.time())
        try:
            # unix only flag disables the fifos and inet connections
            p = Popen([self._ds9_path,
                       "-xpa", "local",
                       "-unix_only", "-title", title,
                       "-unix", "{0:s}".format(unix_name)],
                      shell=False, env=env)

            # wait until ds9 starts and the .IMT socket exists
            while wait_time > 0:
                file_list = os.listdir(self._tmpd_name)
                if ".IMT" in file_list and any(['DS9' in fl for fl in file_list]):
                    break
                time.sleep(0.5)
                wait_time -= 0.5

            if wait_time == 0:
                from signal import SIGTERM
                os.kill(p.pid, SIGTERM)
                print(f"Connection timeout with the ds9. Try to increase the \
                    *wait_time* parameter (current value \
                    is  {self.wait_time} s)")

        except (OSError, ValueError, AttributeError) as e:
            warnings.warn("Starting ds9 failed")
            shutil.rmtree(self._tmpd_name)

        else:
            self._tmp_dir = self._tmpd_name
            self._ds9_process = p
            self._process_list.append(p)

        # this might be sketchy
        try:
            file_list.remove(".IMT")  # should be in the directory, if not
        except (ValueError, IOError):
            warnings.warn("IMT not found in tmp, using first thing in list")
        if len(file_list) > 0:
            for fl in file_list:
                if 'DS9' in fl: break
            xpaname = os.path.join(self._tmpd_name, fl)
        else:
            shutil.rmtree(self._tmpd_name)
            raise ValueError("Problem starting ds9 local socket connection")

        env["XPA_TMPDIR"] = "/tmp/xpa"  # for all local connections
        self._need_to_purge = True
        self._xpa_method = 'local'
        return xpaname, unix_name

    def set_iraf_display(self):
        """Set the environemnt variable IMTDEV to the current display.

        the socket address of the current imexam.ds9 instance is used
        Notes
        -----
        For example, your pyraf commands will use this ds9 for display.

        TODO: Not sure this is still needed. Stop using IRAF.
        """
        os.environ["IMTDEV"] = "unix:{0:s}".format(self._ds9_unix_name)