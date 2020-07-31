import multiprocessing as mp
import time
from collections import Sequence

from luna.util.progress_tracker import ProgressTracker

import logging
logger = logging.getLogger()


MAX_NPROCS = mp.cpu_count() - 1


class Sentinel:
    pass


class ParallelJobs:

    # TODO: add option to Threads/Multiprocessing
    def __init__(self, nproc=MAX_NPROCS):

        self.nproc = nproc

    def _exec_func(self, data, func):
        start = time.time()
        failed = False
        try:
            # Unpack args as a dictionary.
            if isinstance(data, dict):
                func(**data)
            # Unpack args as a list.
            elif isinstance(data, Sequence):
                func(*data)
            # Try unpack the function as a single arg.
            else:
                func(data)
        # Capture any errors while executing the above code.
        except Exception as e:
            logger.exception(e)
            failed = True

        proc_time = time.time() - start

        return proc_time, failed

    def _producer(self, args_list, job_queue):
        for data in args_list:
            job_queue.put(data)

    def _consumer(self, func, job_queue, progress_queue):
        while True:
            data = job_queue.get()

            # If sentinel is found, break.
            if isinstance(data, Sentinel):
                break

            # Execute provided function.
            proc_time, failed = self._exec_func(data, func)
            # Update progress tracker.
            progress_queue.put((data, proc_time, failed))

            job_queue.task_done()

    def _sequential(self, args_list, func, progress_queue):
        # Run jobs sequentially.
        for data in args_list:
            # Execute provided function.
            proc_time, failed = self._exec_func(data, func)
            # Update progress tracker.
            progress_queue.put((data, proc_time, failed))

    def run_jobs(self, args_list, consumer_func, job_name=None):

        # Queue for progress tracker.
        progress_queue = mp.JoinableQueue(maxsize=1)

        # Progress tracker
        job_processing = ProgressTracker(len(args_list), progress_queue, job_name)
        job_processing.start()

        # Initialize a new progress bar (display a 0% progress).
        progress_queue.put(None)

        if self.nproc is not None:
            job_queue = mp.JoinableQueue(maxsize=self.nproc)

            for i in range(self.nproc):
                p = mp.Process(name="ConsumerProcess-%d" % i, target=self._consumer, args=(consumer_func, job_queue, progress_queue))
                p.daemon = True
                p.start()

            # Produce tasks to consumers.
            self._producer(args_list, job_queue)

            # Join all processes and finalize the progress tracker.
            job_queue.join()

            # Finish the progress tracker.
            job_processing.end()

            # Add sentinels to stop consumers.
            sentinel = Sentinel()
            [job_queue.put(sentinel) for i in range(self.nproc)]
        else:
            self._sequential(args_list, consumer_func, progress_queue)

        return job_processing.errors
