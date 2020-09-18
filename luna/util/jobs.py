import multiprocessing as mp
import time
from collections import Sequence

from luna.util.progress_tracker import ProgressData, ProgressTracker
from luna.util.file import get_unique_filename

import logging
logger = logging.getLogger()


MAX_NPROCS = mp.cpu_count() - 1


class Sentinel:
    pass


class ParallelJobs:

    # TODO: add option to Threads/Multiprocessing
    def __init__(self, nproc=MAX_NPROCS):

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

    def _producer(self, args_list, job_queue):
        for data in args_list:
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
            pd = ProgressData(input_data=data, output_data=output, exception=exception, proc_time=proc_time, func=func)
            progress_queue.put(pd)

            job_queue.task_done()

    def _saver(self, output_queue, output_file, proc_func=None, output_header=None):

        with open(output_file, "w") as OUT:
            if output_header is not None:
                OUT.write(output_header.strip())
                OUT.write("\n")

            while True:
                data = output_queue.get()

                # If sentinel is found, break.
                if isinstance(data, Sentinel):
                    break

                line = data
                if proc_func is not None:
                    # Execute the provided function.
                    output, exception, proc_time = self._exec_func(data, proc_func)
                    line = output

                if output is not None:
                    try:
                        OUT.write(str(line).strip())
                        OUT.write("\n")
                        OUT.flush()
                    except Exception as e:
                        logger.error("An error occurred while trying to save the output '%s'." % str(line))
                        logger.exception(e)

                output_queue.task_done()

    def _sequential(self, args_list, func, progress_queue):
        # Run jobs sequentially.
        for data in args_list:
            # Execute provided function.
            output, exception, proc_time = self._exec_func(data, func)

            # Save data.
            pd = ProgressData(input_data=data, output_data=output, exception=exception, proc_time=proc_time)

            # Update progress tracker.
            progress_queue.put(pd)

    def run_jobs(self, args_list, consumer_func, output_file=None, proc_output_func=None, output_header=None, job_name=None):

        if proc_output_func is not None and output_file is None:
            output_file = get_unique_filename(".") + ".output"
            logger.warning("No output file was defined. So, it will try to save results at '%s'." % output_file)
        elif output_file is not None:
            logger.warning("The output file '%s' was defined. So, it will try to save results at it." % output_file)

        # Queue for progress tracker.
        progress_queue = mp.JoinableQueue(maxsize=1)

        # Progress tracker
        self.progress_tracker = ProgressTracker(len(args_list), progress_queue, job_name)
        self.progress_tracker.start()

        # Initialize a new progress bar (display a 0% progress).
        progress_queue.put(None)

        if self.nproc is not None:
            job_queue = mp.JoinableQueue(maxsize=self.nproc)

            output_queue = None
            if output_file is not None:
                output_queue = mp.JoinableQueue()

            for i in range(self.nproc):
                p = mp.Process(name="ConsumerProcess-%d" % i, target=self._consumer, args=(consumer_func, job_queue, progress_queue,
                                                                                           output_queue,))
                p.daemon = True
                p.start()

            if output_file is not None:
                o = mp.Process(name="WriterProcess-%d" % i, target=self._saver, args=(output_queue, output_file,
                                                                                      proc_output_func, output_header,))
                o.daemon = True
                o.start()

            # Produce tasks to consumers.
            self._producer(args_list, job_queue)

            # Sentinels to stop consumers.
            sentinel = Sentinel()

            # Join all processes and add sentinels to stop consumers.
            job_queue.join()
            [job_queue.put(sentinel) for i in range(self.nproc)]

            if output_file is not None:
                output_queue.join()
                output_queue.put(sentinel)

        else:
            self._sequential(args_list, consumer_func, progress_queue)

        # Finish the progress tracker.
        self.progress_tracker.end()

        return self.progress_tracker.results
