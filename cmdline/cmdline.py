#!/usr/bin/env python

"""Base class for building command line applications.

This class includes logging, exception handling, and argument processing.  It
also sets up an option parser for the user to simply extend.

Version 1.1
Version 1.0 was by ???

Copyright (c) 2009, John Wiegley.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

- Neither the name of New Artisans LLC nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import sys
import os
import time
import optparse
import inspect
import logging
import logging.handlers

LEVELS = {'DEBUG':    logging.DEBUG,
          'INFO':     logging.INFO,
          'WARNING':  logging.WARNING,
          'ERROR':    logging.ERROR,
          'CRITICAL': logging.CRITICAL}

class CommandLineApp(object):
    "Base class for building command line applications."

    force_exit  = True           # If true, always ends run() with sys.exit()
    log_handler = None

    options = {
        'debug':    False,
        'verbose':  False,
        'logfile':  False,
        'loglevel': False
    }

    def __init__(self):
        "Initialize CommandLineApp."
        # Create the logger
        self.log = logging.getLogger(os.path.basename(sys.argv[0]))
        ch = logging.StreamHandler()
        formatter = logging.Formatter("%(name)s: %(levelname)s: %(message)s")
        ch.setFormatter(formatter)
        self.log.addHandler(ch)
        self.log_handler = ch

        # Setup the options parser
        usage = 'usage: %prog [options] <BOUND-IP-ADDRESS>'
        op = self.option_parser = optparse.OptionParser(usage = usage)

        op.add_option('', '--debug',
                      action='store_true', dest='debug',
                      default=False, help='show debug messages and pass exceptions')
        op.add_option('-v', '--verbose',
                      action='store_true', dest='verbose',
                      default=False, help='show informational messages')
        op.add_option('-q', '--quiet',
                      action='store_true', dest='quiet',
                      default=False, help='do not show log messages on console')
        op.add_option('', '--log', metavar='FILE',
                      type='string', action='store', dest='logfile',
                      default=False, help='append logging data to FILE')
        op.add_option('', '--loglevel', metavar='LEVEL',
                      type='string', action='store', dest='loglevel',
                      default=False, help='set log level: DEBUG, INFO, WARNING, ERROR, CRITICAL')
        return

    def main(self, *args):
        """Main body of your application.

        This is the main portion of the app, and is run after all of the
        arguments are processed.  Override this method to implment the primary
        processing section of your application."""
        pass

    def handleInterrupt(self):
        """Called when the program is interrupted via Control-C or SIGINT.
        Returns exit code."""
        self.log.error('Canceled by user.')
        return 1

    def handleMainException(self):
        "Invoked when there is an error in the main() method."
        if not self.options.debug:
            self.log.exception('Caught exception')
        return 1

    ## INTERNALS (Subclasses should not need to override these methods)

    def run(self):
        """Entry point.

        Process options and execute callback functions as needed.  This method
        should not need to be overridden, if the main() method is defined."""
        # Process the options supported and given
        self.options, main_args = self.option_parser.parse_args()

        if self.options.logfile:
            fh = logging.handlers.RotatingFileHandler(self.options.logfile,
                                                      maxBytes = (1024 * 1024),
                                                      backupCount = 5)
            formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
            fh.setFormatter(formatter)
            self.log.addHandler(fh)

        if self.options.quiet:
            self.log.removeHandler(self.log_handler)
            ch = logging.handlers.SysLogHandler()
            formatter = logging.Formatter("%(name)s: %(levelname)s: %(message)s")
            ch.setFormatter(formatter)
            self.log.addHandler(ch)
            self.log_handler = ch

        if self.options.loglevel:
            self.log.setLevel(LEVELS[self.options.loglevel])
        elif self.options.debug:
            self.log.setLevel(logging.DEBUG)
        elif self.options.verbose:
            self.log.setLevel(logging.INFO)
        
        exit_code = 0
        try:
            # We could just call main() and catch a TypeError, but that would
            # not let us differentiate between application errors and a case
            # where the user has not passed us enough arguments.  So, we check
            # the argument count ourself.
            argspec = inspect.getargspec(self.main)
            expected_arg_count = len(argspec[0]) - 1

            if len(main_args) >= expected_arg_count:
                exit_code = self.main(*main_args)
            else:
                self.log.debug('Incorrect argument count (expected %d, got %d)' %
                               (expected_arg_count, len(main_args)))
                self.option_parser.print_help()
                exit_code = 1

        except KeyboardInterrupt:
            exit_code = self.handleInterrupt()

        except SystemExit, msg:
            exit_code = msg.args[0]

        except Exception:
            exit_code = self.handleMainException()
            if self.options.debug:
                raise
            
        if self.force_exit:
            sys.exit(exit_code)
        return exit_code

    def http_uptest(self, proc, hostname = 'localhost', port = 80, url = u'/'):
        import httplib
        conn = httplib.HTTPConnection(hostname, port = port)
        try:
            conn.request("GET", url)
        except Exception:
            self.log.exception('-- http_uptest exception:')
            return False

        resp = conn.getresponse()
        conn.close()
        if resp.status == 200:
            self.log.debug('-- http_uptest succeeded: %s' %
                           (resp.status, resp.reason))
            return True
        else:
            self.log.warning('-- http_uptest FAILED: %s %s' %
                             (resp.status, resp.reason))
            return False

    def spawn_and_wait(self, cmd, *args, **kwargs):
        from subprocess import Popen

        p = Popen((cmd,) + args)

        while True:
            # Wait a bit before checking on the process
            if kwargs.has_key('poll'):
                time.sleep(kwargs['poll'])
            else:
                time.sleep(1)

            # Check whether the process aborted entirely for any reason.  If
            # so, log the fact and then let our outer loop run it again.
            sts = p.poll()
            if sts is not None:
                self.log.info('-- %s exited: %d' % (cmd, sts))
                return sts

            # The process is still running.  Check whether it is still viable
            # by calling the given callback.
            death = False
            try:
                if kwargs.has_key('uptest') and \
                   callable(kwargs['uptest']) and \
                   not kwargs['uptest'](p):
                    death = True

            except Exception:
                self.log.exception('-- %s exception:' % cmd)

            # If the process is no longer viable, we kill it and exit
            if death is True:
                try:
                    import win32api
                    import win32con
                    import win32process

                    handle = win32api.OpenProcess(win32con.PROCESS_ALL_ACCESS,
                                                  True, p.pid)
                    exitcode = win32process.GetExitCodeProcess(handle)
                    if exitcode == win32con.STILL_ACTIVE:
                        win32api.TerminateProcess(handle, 0)
                        self.log.warning('-- %s killed' % cmd)
                except:
                    import signal
                    try: os.kill(p.pid, signal.SIGHUP)
                    except: pass
                    try: os.kill(p.pid, signal.SIGINT)
                    except: pass
                    try: os.kill(p.pid, signal.SIGQUIT)
                    except: pass
                    try: os.kill(p.pid, signal.SIGKILL)
                    except: pass
                    self.log.warning('-- %s killed' % cmd)

                return -1


if __name__ == "__main__":
    import tempfile

    class sample(CommandLineApp):
        def main(self, *args):
            self.log.info("Saw args: %s" % (args,))
            temp = tempfile.NamedTemporaryFile()
            temp.write("Args were: %s" % (args,))
            temp.flush()
            name = temp.name
            self.log.info("Wrote args to '%s'" % name)
            self.spawn_and_wait('/bin/ls', '-l', name)
            del temp
            self.spawn_and_wait('/bin/ls', '-l', name)
            self.spawn_and_wait('/bin/sleep', '30', uptest = self.http_uptest)
            x = 1 / 0           # logged as an exception

    sample().run()
