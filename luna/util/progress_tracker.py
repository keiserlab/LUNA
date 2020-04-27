# Source: http://www.codekoala.com/posts/command-line-progress-bar-python/

from queue import Queue
from threading import Thread, Event

import time
import sys


import logging
logger = logging.getLogger()


class ProgressTracker:

    def __init__(self, ntasks, task_name=None):
        self.ntasks = ntasks
        self.task_name = task_name

        self._progress = 0

        self.queue = Queue(1)   # used to communicate progress to the thread
        self.event = Event()    # used to tell the thread when to finish
        self.progress_bar = Thread(target=self.print_progress, args=(self.event, self.queue))
        self.progress_bar.daemon = True

        # Save any errors found during the task processing.
        self.errors = []
        self.running_times = []
        self._start_time = None
        self._end_time = None

    def print_progress(self, e, q):
        """Updates a progress bar on stdout anytime progress is made"""

        while True:
            while not q.full():
                # Stops the loop if function end() is called before the queue gets full.
                if e.is_set():
                    break
                # wait for more progress to be made
                time.sleep(0.1)

            # If our event is set and there is any progress in the queue,
            # break out of the infinite loop and prepare to terminate this thread
            if e.is_set() and q.full() is False:
                break

            # get the current progress value
            p = q.get()
            perc = round((p / self.ntasks), 2) * 100

            task_name = ""
            if self.task_name:
                task_name = " - %s" % self.task_name

            # print the current progress value in bar and percent form
            sys.stdout.write('\r%s%% [%s] %d/%d [Avg: %.2fs/task]%s' % (int(perc), ("\u25A0" * int(perc / 2)).ljust(50, ' '), p,
                                                                        self.ntasks, self.avg_running_time, task_name))
            sys.stdout.flush()

    @property
    def progress(self):
        """Returns the current progress value"""
        return self._progress

    @progress.setter
    def progress(self, value):
        """Sets the current progress value, passing updates to the thread"""
        self._progress = value

        # this is the key to this solution
        self.queue.put(self._progress)

    @property
    def running_time(self):
        if self._start_time and self._end_time:
            return round(self._end_time - self._start_time, 2)
        return None

    @property
    def avg_running_time(self):
        return round(sum(self.running_times) / len(self.running_times), 2)

    def start(self):
        self._start_time = round(time.time(), 2)
        self.progress_bar.start()

    def end(self):
        self.event.set()
        self.progress_bar.join()

        self._end_time = round(time.time(), 2)

        sys.stdout.write('\n\n')
