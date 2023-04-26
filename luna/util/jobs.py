import multiprocessing as mp
import time
from collections import Sequence

from luna.util.progress_tracker import ProgressData, ProgressTracker
from luna.util.file import new_unique_filename

import logging
logger = logging.getLogger()


MAX_NPROCS = mp.cpu_count() - 1


class Sentinel:
    """Custom sentinel to stop workers"""
    pass


class ArgsGenerator:
    """Custom generator that implements __len__().
       This class can be used in conjunction with
       :class:`~luna.util.progress_tracker.ProgressTracker` in cases
       where the tasks are obtained from generators. Note that
       :class:`~luna.util.progress_tracker.ProgressTracker` requires a
       pre-defined number of tasks to calculate the progress, therefore a
       standard generator cannot be used directly as it does not implement
       __len__(). Then, with `ArgsGenerator`, one may take advantage of
       generators and :class:`~luna.util.progress_tracker.ProgressTracker`
       by explicitly providing the number of tasks that will be generated.

       Parameters
       ----------
       generator : generator
           The tasks generator.
       nargs : int
           The number of tasks that will be generated.
    """

    def __init__(self, generator, nargs):
        self.generator = generator
        self.nargs = nargs

    def __len__(self):
        return self.nargs

    def __iter__(self):
        for d in self.generator:
            yield d


class ParallelJobs:

    """Executes a set of tasks in parallel
    (:py:class:`~multiprocessing.JoinableQueue`) or sequentially.

    Parameters
    ----------
    nproc : int or None
       The number of CPUs to use. The default value is the ``maximum number
       of CPUs - 1``. If ``nproc`` is None, 0, or 1, run the jobs sequentially.
       Otherwise, use the ``maximum number of CPUs - 1``.

    Attributes
    ----------
    nproc : int
        The number of CPUs to use.
    progress_tracker : ProgressTracker
        A :class:`~luna.util.progress_tracker.ProgressTracker` object to track
        the tasks' progress.
    """

    # TODO: add option to Threads/Multiprocessing
    def __init__(self, nproc=MAX_NPROCS):

        if nproc is not None:
            # Use 'MAX_NPROCS' if a non-integer has been provided.
            if not isinstance(nproc, int) or nproc < 0:
                nproc = MAX_NPROCS
            elif nproc in [0, 1]:
                nproc = None
            else:
                nproc = min(nproc, MAX_NPROCS)

        self.nproc = nproc
        self.progress_tracker = None

    def _exec_func(self, data, func):
        start = time.time()

        output = None
        exception = None
        try:
            # Unpack args as a dictionary.
            if isinstance(data, dict):
                output = func(**data)
            # Unpack args as a list.
            elif isinstance(data, Sequence):
                output = func(*data)
            # Try unpack the function as a single arg.
            else:
                output = func(data)
        # Capture any errors while executing the above code.
        except Exception as e:
            logger.exception(e)
            exception = e

        proc_time = time.time() - start

        return output, exception, proc_time

    def _producer(self, args, job_queue):
        for data in args:
            job_queue.put(data)

    def _consumer(self, func, job_queue, progress_queue, output_queue=None):
        while True:
            data = job_queue.get()

            # If sentinel is found, break.
            if isinstance(data, Sentinel):
                break

            # Execute the provided function.
            output, exception, proc_time = self._exec_func(data, func)

            if output is not None and output_queue is not None:
                output_queue.put((data, output))

            # Update progress tracker.
            pd = ProgressData(input_data=data,
                              output_data=output,
                              exception=exception,
                              proc_time=proc_time,
                              func=func)
            progress_queue.put(pd)

            job_queue.task_done()

    def _saver(self, output_queue, output_file,
               proc_func=None, output_header=None):

        with open(output_file, "w") as OUT:
            if output_header is not None:
                OUT.write(output_header.strip())
                OUT.write("\n")

            while True:
                data = output_queue.get()

                # If sentinel is found, break.
                if isinstance(data, Sentinel):
                    break

                line = None
                if proc_func is not None:
                    # Execute the provided function.
                    output, exception, proc_time = self._exec_func(data,
                                                                   proc_func)
                    line = output

                try:
                    # If no data is stored in line, try to access the
                    # output generated by the _consumer() function.
                    if line is None:
                        line = data[1]

                    OUT.write(str(line).strip())
                    OUT.write("\n")
                    OUT.flush()
                except Exception as e:
                    logger.error("An error occurred while trying to save "
                                 "the output '%s'." % str(line))
                    logger.exception(e)

                output_queue.task_done()

    def _sequential(self, args, func, progress_queue):
        # Run jobs sequentially.
        for data in args:
            # Execute provided function.
            output, exception, proc_time = self._exec_func(data, func)

            # Save data.
            pd = ProgressData(input_data=data,
                              output_data=output,
                              exception=exception,
                              proc_time=proc_time)

            # Update progress tracker.
            progress_queue.put(pd)

    def run_jobs(self, args, consumer_func, output_file=None,
                 proc_output_func=None, output_header=None, job_name=None):
        """
        Run a set of tasks in parallel or sequentially according to the
        ``nproc``.

        Parameters
        ----------
        args : iterable of iterables, `ArgsGenerator`
            A sequence of arguments to be provided to the consumer function
            ``consumer_func``.
        consumer_func : function
            The function that will be executed for each set of arguments in
            ``args``.
        output_file : str, optional
            Save outputs to this file.
            If ``proc_output_func`` is not provided, it tries to save a
            stringified version of each output data. Otherwise, it executes
            ``proc_output_func`` first and its output will be printed to the
            output file instead.

            Note: if ``proc_output_func`` is provided but not ``output_file``,
            a new random unique filename will be generated and the file will be
            saved in the current directory.
        proc_output_func : function, optional
            Post-processing function that is executed for each output data
            produced by ``consumer_func``.
        output_header : str, optional
            A header for the output file.
        job_name : str, optional
            A name to identify the job.

        Returns
        -------
         : :class:`~luna.util.progress_tracker.ProgressResult`
        """
        if proc_output_func is not None and output_file is None:
            output_file = new_unique_filename(".") + ".output"
            logger.warning("No output file was defined. So, it will try to "
                           "save results at '%s'." % output_file)
        elif output_file is not None:
            logger.warning("The output file '%s' was defined. So, it will "
                           "try to save results at it." % output_file)

        # Queue for progress tracker.
        progress_queue = mp.JoinableQueue(maxsize=1)

        # Progress tracker
        self.progress_tracker = ProgressTracker(len(args), progress_queue,
                                                job_name)
        self.progress_tracker.start()

        # Initialize a new progress bar (display a 0% progress).
        progress_queue.put(None)

        if self.nproc is not None:
            job_queue = mp.JoinableQueue(maxsize=self.nproc)

            output_queue = None
            if output_file is not None:
                output_queue = mp.JoinableQueue()

            for i in range(self.nproc):
                p = mp.Process(name="ConsumerProcess-%d" % i,
                               target=self._consumer,
                               args=(consumer_func, job_queue,
                                     progress_queue, output_queue,))
                p.daemon = True
                p.start()

            if output_file is not None:
                o = mp.Process(name="WriterProcess-%d" % i,
                               target=self._saver,
                               args=(output_queue, output_file,
                                     proc_output_func, output_header,))
                o.daemon = True
                o.start()

            # Produce tasks to consumers.
            self._producer(args, job_queue)

            # Sentinels to stop consumers.
            sentinel = Sentinel()

            # Join all processes and add sentinels to stop consumers.
            job_queue.join()
            [job_queue.put(sentinel) for i in range(self.nproc)]

            if output_file is not None:
                output_queue.join()
                output_queue.put(sentinel)

        else:
            self._sequential(args, consumer_func, progress_queue)

        # Finish the progress tracker.
        self.progress_tracker.end()

        return self.progress_tracker.results
