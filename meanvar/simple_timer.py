import time
import sys

class SimpleTimer:
    """Simple abstraction to allow for basic timing. """

    def __init__(self):
        self.start = time.time()

    def report(self, msg, do_reset=False, file=sys.stdout):
        """Print to stdout msg followed by the runtime.

        When true, do_reset will result in a reset of start time.

        """
        print >> file, "%s (%s s)" % (msg, time.time() - self.start)
        if do_reset:
            self.start = time.time()

    def result(self, msg, do_reset=False):
        """Return log message containing ellapsed time as a string.

        When true, do_reset will result in a reset of start time.

        """
        result = "%s (%s s)" % (msg, time.time() - self.start)
        if do_reset:
            self.start = time.time()
        return result

    def reset(self):
        """Reset start time"""
        self.start = time.time()


    def runtime(self):
        """Return ellapsed time and reset start. """
        t = time.time() - self.start
        self.start = time.time()
        return t
