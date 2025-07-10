"""
Run commands on multiple processors in python.

Adapted from examples on https://docs.python.org/2/library/multiprocessing.html
"""

# These functions are also borrow from the Gypsum-DL script Parallelizer.py
# These functions were renamed to be pep8 compliant
# ie )
# def multi_threading became def multi_threading


import multiprocessing
from typing import Callable


def multi_threading(inputs: list, num_processors: int, task_name: Callable) -> list:
    """Initialize this object.

    Args:
        inputs ([data]): A list of data. Each datum contains the details to
            run a single job on a single processor.
        num_processors (int): The number of processors to use.
        task_name (class): The class that governs what to do for each
            job on each processor.
    """

    results = []

    # If there are no inputs, just return an empty list.
    if not inputs:
        return results

    num_processors = count_processors(len(inputs), num_processors)

    tasks = []

    for index, item in enumerate(inputs):
        if not isinstance(item, tuple):
            item = (item,)
        task = (index, (task_name, item))
        tasks.append(task)

    if num_processors == 1:
        for item in tasks:
            job, args = item[1]
            output = job(*args)
            results.append(output)
    else:
        results = start_processes(tasks, num_processors)

    return results


# Worker function


def worker(input: multiprocessing.Queue, output: multiprocessing.Queue) -> None:
    """
    Worker function to process tasks in parallel.

    Args:
        input (multiprocessing.Queue): Queue containing tasks to execute.
        output (multiprocessing.Queue): Queue to store task results.
    """
    for seq, job in iter(input.get, "STOP"):
        func, args = job
        result = func(*args)
        ret_val = (seq, result)
        output.put(ret_val)


def count_processors(num_inputs: int, num_processors: int) -> int:
    """
    Checks processors available and returns a safe number of them to
    utilize.

    :param int num_inputs: The number of inputs.
    :param int num_processors: The number of desired processors.

    :returns: The number of processors to use.
    """
    # first, if num_processors<= 0, determine the number of processors to
    # use programatically
    if num_processors <= 0:
        num_processors = multiprocessing.cpu_count()

    # reduce the number of processors if too many have been specified
    if num_inputs < num_processors:
        num_processors = num_inputs

    return num_processors


def start_processes(inputs: list, num_processors: int) -> list:
    """
    Creates and starts worker processes to handle tasks in parallel.

    Args:
        inputs (list): A list of tasks to be processed.
        num_processors (int): The number of processors to use.

    Returns:
        list: The results of the tasks processed in parallel.
    """
    # Create queues
    task_queue = multiprocessing.Queue()
    done_queue = multiprocessing.Queue()

    # Submit tasks
    for item in inputs:
        task_queue.put(item)

    # Start worker processes
    for _ in range(num_processors):
        multiprocessing.Process(target=worker, args=(task_queue, done_queue)).start()

    results = [done_queue.get() for _ in range(len(inputs))]
    # Tell child processes to stop
    for _ in range(num_processors):
        task_queue.put("STOP")

    results.sort(key=lambda tup: tup[0])

    return [item[1] for item in map(list, results)]


###
# Helper functions
###


def flatten_list(tier_list: list[list]) -> list:
    """
    Given a list of lists, this returns a flat list of all items.

    :params list tier_list: A 2D list.

    :returns: A flat list of all items.
    """
    if tier_list is None:
        return []
    return [item for sublist in tier_list for item in sublist]


def strip_none(none_list: list) -> list:
    """
    Given a list that might contain None items, this returns a list with no
    None items.

    :params list none_list: A list that may contain None items.

    :returns: A list stripped of None items.
    """
    return [] if none_list is None else [x for x in none_list if x is not None]
