# Source: http://www.codekoala.com/posts/command-line-progress-bar-python/

from threading import Thread, Event

import time
import sys
from colorlog.escape_codes import escape_codes, parse_colors


import logging
logger = logging.getLogger()


class ProgressData:
    """
    Custom data structure to store data generated during the
    execution of a single task.

    Parameters
    ----------
    input_data : any
        The input data of the executed task.
    proc_time : float
        How long a task took to be executed.
    output_data : any, optional
        The output data produced by a given task.
        The data type is the same as the executed function's return.
    exception : :py:class:`Exception`, optional
        If an exception was raised, then ``exception`` stores an
        Exception object. Otherwise, ``exception`` will be set to None.
    func : function, optional
        The executed function for reference.
    """

    def __init__(self, input_data, proc_time,
                 output_data=None, exception=None, func=None):

        self.input_data = input_data
        self.output_data = output_data
        self.proc_time = proc_time
        self.exception = exception
        self.func = func


class ProgressResult:
    """
    Custom iterable class to store `ProgressData` objects as they are produced
    during the execution of a set of tasks.

    This class implements append() and __iter__() by default.

    Parameters
    ----------
    results : list, optional
        A pre-populated list of `ProgressData` objects.

    Attributes
    ----------
    results : list
        The list of `ProgressData` objects.
    """

    def __init__(self, results=None):
        if results and not isinstance(results, list):
            raise TypeError("A list was expected but a different data type "
                            "was provided.")
        self.results = results or []

    @property
    def inputs(self):
        """list: Inputs from each `ProgressData` object stored in \
        ``results``."""
        try:
            return [r.input_data for r in self.results]
        except AttributeError:
            raise TypeError("An invalid object was found in 'results'. "
                            "Only ProgressData objects are valid.")

    @property
    def outputs(self):
        """list of tuple: Outputs from each `ProgressData` object stored in
        ``results``. Each tuple contains the input and the output produced for
        that input."""
        try:
            return [(r.input_data, r.output_data) for r in self.results]
        except AttributeError:
            raise TypeError("An invalid object was found in 'results'. "
                            "Only ProgressData objects are valid.")

    @property
    def errors(self):
        """list of tuple: Errors from each `ProgressData` object stored in \
        ``results``.
        Each tuple contains the input and the exception raised during the
        execution of a task with that input."""
        try:
            return [(r.input_data, r.exception) for r in self.results
                    if r.exception is not None]
        except AttributeError:
            raise TypeError("An invalid object was found in 'results'. "
                            "Only ProgressData objects are valid.")

    def append(self, r):
        """Add a new `ProgressData` object to ``results``"""
        if not isinstance(r, ProgressData):
            raise TypeError("Only ProgressData objects are valid.")
        self.results.append(r)

    def __len__(self):
        return len(self.results)

    def __iter__(self):
        for r in self.results:
            yield r


class ProgressTracker:
    """
    A progress tracker for tasks queued in ``queue``.

    Parameters
    ----------
    ntasks : int
        The number of tasks to be executed.
    queue : :py:class:`~multiprocessing.Queue`
        A queue to track the tasks' progress.
    task_name : str, optional
        A name to identify the set of tasks.

    Attributes
    ----------
    ntasks : int
        The number of tasks to be executed.
    queue : :py:class:`~multiprocessing.Queue`
        The queue to track the tasks' progress from.
    task_name : str, optional
        The name to identify the set of tasks.
    results : ProgressResult
        Store results and any errors found during the tasks' processing.
    nerrors : int
        Number of errors.
    """

    def __init__(self, ntasks, queue, task_name=None):
        self.ntasks = ntasks
        self.task_name = task_name

        self._progress = 0

        self.queue = queue   # used to communicate progress to the thread
        self._event = Event()    # used to tell the thread when to finish
        self._progress_bar = Thread(target=self._print_progress,
                                    args=(self._event, self.queue))
        self._progress_bar.daemon = True

        # Save results and any errors found during the task processing.
        self.results = ProgressResult()
        self.nerrors = 0

        self._running_times = []
        self._start_time = None
        self._end_time = None

    def _show_progress_bar(self, p, perc):
        task_name = ""
        if self.task_name:
            task_name = " - %s" % self.task_name

        msg = ('%s%% [%s] %d/%d [Avg: %.2fs/task; Errors: %d]%s.'
               % (int(perc), ("\u25A0" * int(perc / 2)).ljust(50, ' '),
                  p, self.ntasks, self.avg_running_time,
                  self.nerrors, task_name))

        format_str = '\r[%s]    %s%s %s%s  %s'
        progress_str = \
            (format_str % (time.strftime('%Y-%m-%d %H:%M:%S'),
             parse_colors("purple"),
             "PROGRESS".ljust(10, " "),
             escape_codes["reset"],
             "".rjust(26, " "),
             msg))

        sys.stdout.write(progress_str)
        sys.stdout.flush()

    def _print_progress(self, e, q):
        """Updates a progress bar on stdout anytime progress is made"""

        while True:
            while not q.full():
                # Stops the loop if function end() is called before the queue
                # gets full.
                if e.is_set():
                    break
                # wait for more progress to be made
                time.sleep(0.1)

            # If our event is set and there is any progress in the queue,
            # break out of the infinite loop and prepare to terminate this
            # thread
            if e.is_set() and q.full() is False:
                break

            # get the current progress data.
            progress_data = q.get()

            if progress_data is not None:
                self.results.append(progress_data)

                self._running_times.append(progress_data.proc_time)

                if progress_data.exception is not None:
                    self.nerrors += 1

                self.progress += 1

            perc = round((self.progress / self.ntasks), 2) * 100 if self.ntasks > 0 else 0
            self._show_progress_bar(self.progress, perc)

    @property
    def progress(self):
        """int: Current number of tasks executed."""
        return self._progress

    @progress.setter
    def progress(self, value):
        self._progress = value

    @property
    def running_time(self):
        """float: Total running time."""
        if self._start_time and self._end_time:
            return round(self._end_time - self._start_time, 2)
        return None

    @property
    def avg_running_time(self):
        """float: Average running time."""
        if self._running_times:
            total = sum(self._running_times)
            count = len(self._running_times)
            return round(total / count, 2)
        else:
            return 0

    def start(self):
        """ Start the progress tracker. """
        self._start_time = round(time.time(), 2)
        self._progress_bar.start()

    def end(self):
        """ Finish the progress tracker. """
        self._event.set()
        self._progress_bar.join()
        self._end_time = round(time.time(), 2)
        sys.stdout.write('\n')
