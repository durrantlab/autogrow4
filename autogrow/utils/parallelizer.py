# Copyright 2018 Jacob D. Durrant

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# NOTE: This is for parallelizing Python code. For command-line commands, use a
# plugin. Adapted from Gypsum-DL codebase:
# https://github.com/durrantlab/gypsum_dl/blob/1.2.1/gypsum_dl/Parallelizer.py

"""
Abstract parallel computation utility.

The "parallelizer" object exposes a simple map interface that takes a function
and a list of arguments and returns the result of applying the function to
each argument. Internally, the parallelizer class can determine what parallel
capabilities are present on a system and automatically pick between
"multiprocessing" or "serial" in order to speed up the map operation. This
approach simplifies development and allows the same program to run on a laptop
or a high-performance computer cluster, utilizing the full resources of each
system. (Description provided by Harrison Green.)
"""
import __future__
import multiprocessing


class Parallelizer(object):
    """Abstract class for parallel computation management.

    Exposes a simple map interface that takes a function and list of arguments
    and returns the result of applying the function to each argument.
    Automatically detects system parallel capabilities and chooses between
    'multiprocessing' or 'serial' modes to optimize performance.

    Attributes:
        mode (str): Mode used for parallelization. Determined by mode parameter
            and environment. Defaults to multiprocessing unless specified
            otherwise.
        parallel_obj: Instantiated parallelization class object. Set to None for
            simpler methods like serial.
        num_procs (int): Number of processors/nodes to use. If None, uses all
            available. Forced to 1 if mode is serial.
    """

    def __init__(self, mode=None, num_procs=None):
        """
        Initialize parallelizer and configures specific multiprocessing class.

        Args:
            mode (str, optional): Multiprocess mode to use ('serial',
                'multiprocessing', or None). If None, attempts to choose an
                available multiprocessing option. Best practice is to specify mode
                explicitly, especially for sub-programs used by larger programs.
                Defaults to None.
            num_procs (int, optional): Number of processors/nodes to use. If None,
                uses all available. Overridden to 1 if mode is serial. Defaults to
                None.
        """
        if mode in ["none", "None"]:
            mode = None

        if (
            mode is None
            or mode == "multiprocessing"
            or mode not in ["Serial", "serial"]
        ):
            self.mode = "multiprocessing"

        else:
            self.mode = "serial"

        self.parallel_obj = None

        if self.mode == "serial":
            self.num_procs = 1

        elif num_procs is None or num_procs < 1:
            self.num_procs = self._compute_nodes()

        else:
            self.num_procs = num_procs

    def run(self, args, func, num_procs=None, mode=None):
        """
        Run a task in parallel across the system.

        Args:
            args (list): List of tuples/lists, each containing all parameters
                required by func for a single task
            func (callable): Function object to be executed
            num_procs (int, optional): Number of processors to use. If None, uses
                value from initialization. Cannot override in serial mode. Defaults
                to None.
            mode (str, optional): Multiprocess mode to use. Must be
                'multiprocessing' or 'serial'. If None, uses mode from
                initialization. Primarily for developers. Best to leave as None.
                Defaults to None.

        Returns:
            list: Results from all parallel processes

        Raises:
            Exception: If incompatible mode specified or attempting to override
                num_procs in serial mode

        Example:
            To multiprocess function foo(x,y) for all permutations of x,y from 0-2:
            >>> args = [(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)]
            >>> func = foo  # The namespace of foo
            >>> parallelizer.run(args, func)

        Note:
            For sub-programs used by larger programs with multiprocessing,
            explicitly specify mode to avoid conflicts.
        """
        # determine the mode
        if mode is None:
            mode = self.mode
        elif self.mode != mode:
            if mode not in ["serial", "multiprocessing"]:
                printout = (
                    "Overriding function with a multiprocess mode which doesn't match: "
                    + mode
                )
                raise Exception(printout)

        if num_procs is None:
            num_procs = self.num_procs

        if num_procs != self.num_procs and mode == "serial":
            printout = "Can't override num_procs in serial mode"
            raise Exception(printout)

        if mode != self.mode:
            printout = (
                f"changing mode from {mode} to {self.mode} for development purpose"
            )
            print(printout)

        # compute
        if mode == "multiprocessing":
            return _multi_threading(args, num_procs, func)
        else:
            # serial is running the ParallelThreading with num_procs=1
            return _multi_threading(args, 1, func)

    def _compute_nodes(self, mode=None):
        """
        Compute available compute nodes for selected mode.

        For multiprocessing, returns number of available cores. For serial mode,
        returns 1.

        Args:
            mode (str, optional): Mode to check ('multiprocessing' or 'serial').
                If None, uses current mode. Defaults to None.

        Returns:
            int: Number of compute nodes/processors available
        """
        if mode is None:
            mode = self.mode

        return multiprocessing.cpu_count() if mode == "multiprocessing" else 1

    def return_node(self):
        """
        Return the number of "compute nodes" according to the selected mode.

        For multiprocessing this is the number of available cores
        For serial, this value is 1
        Returns:
        :returns: int num_procs: the number of nodes/processors which is to be used
        """
        return self.num_procs


def _multi_threading(inputs, num_procs, task_name):
    """
    Execute parallel processing using multiprocessing.

    If num_procs is 1, executes tasks serially. Otherwise distributes tasks
    across specified number of processors.

    Args:
        inputs (list): List of data where each item contains details for a
            single job on a single processor
        num_procs (int): Number of processors to use task_name (callable):
        Function that governs processing for each job

    Returns:
        list: Results from all processors. Empty list if no inputs provided.

    Note:
        Automatically formats inputs to tuples and adjusts processor count if
        needed.
    """
    results = []

    # If there are no inputs, just return an empty list.
    if len(inputs) == 0:
        return results

    inputs = _check_and_format_inputs_to_list_of_tuples(inputs)

    num_procs = _count_processors(len(inputs), num_procs)

    tasks = []

    for index, item in enumerate(inputs):
        if not isinstance(item, tuple):
            item = (item,)
        task = (index, (task_name, item))
        tasks.append(task)

    if num_procs == 1:
        for item in tasks:
            job, args = item[1]
            output = job(*args)
            results.append(output)
    else:
        results = _start_processes(tasks, num_procs)

    return results


def _worker(input_q: multiprocessing.Queue, output_q: multiprocessing.Queue):
    """
    Worker function for parallel processing.

    Continuously processes jobs from input queue until receiving 'STOP' signal.

    Args:
        input_q (multiprocessing.Queue): Queue containing (sequence, job)
            pairs. Jobs consist of (function, arguments) pairs.
        output_q (multiprocessing.Queue): Queue for storing (sequence, result)
            pairs.
    """
    for seq, job in iter(input_q.get, "STOP"):
        func, args = job
        result = func(*args)
        ret_val = (seq, result)
        output_q.put(ret_val)

def _check_and_format_inputs_to_list_of_tuples(args):
    """
    Validate and format input arguments list.

    Args:
        args (list): List of input arguments. All items must be the same type
            (either lists or tuples).

    Returns:
        list: List of argument tuples (converts lists to tuples if needed)

    Raises:
        Exception: If args is not a list/tuple, if items are not all same type,
            or if items are neither lists nor tuples
    """
    # Make sure args is a list of tuples
    if type(args) not in [list, tuple]:
        printout = "args must be a list of tuples"
        print(printout)
        raise Exception(printout)

    item_type = type(args[0])
    for i in range(len(args)):
        if type(args[i]) == item_type:
            continue
        printout = "all items within args must be the same type and must be either a list or tuple"
        print(printout)
        raise Exception(printout)
    if item_type == tuple:
        return args
    elif item_type == list:
        args = [tuple(x) for x in args]
        return args
    else:
        printout = "all items within args must be either a list or tuple"
        print(printout)
        raise Exception(printout)


def _count_processors(num_inputs, num_procs):
    """
    Determine appropriate number of processors to use.

    If num_procs <= 0, automatically determines count using multiprocessing.
    Reduces processor count if it exceeds number of inputs.

    Args:
        num_inputs (int): Number of inputs to process
        num_procs (int): Desired number of processors

    Returns:
        int: Number of processors to use
    """
    # first, if num_procs <= 0, determine the number of processors to
    # use programatically
    if num_procs <= 0:
        num_procs = multiprocessing.cpu_count()

    # reduce the number of processors if too many have been specified
    if num_inputs < num_procs:
        num_procs = num_inputs

    return num_procs


def _start_processes(inputs, num_procs):
    """
    Create and manage multiprocessing queues and worker processes.

    Creates input and output queues, starts worker processes, collects results,
    and sends stop signals to workers.

    Args:
        inputs (list): List of (index, (task_name, arguments)) pairs
        num_procs (int): Number of worker processes to create

    Returns:
        list: Results from all workers, sorted by original input order
    """
    # Create queues
    task_queue = multiprocessing.Queue()
    done_queue = multiprocessing.Queue()

    # Submit tasks
    for item in inputs:
        task_queue.put(item)

    # Start worker processes
    for _ in range(num_procs):
        multiprocessing.Process(target=_worker, args=(task_queue, done_queue)).start()
    results = [done_queue.get() for _ in range(len(inputs))]
    # Tell child processes to stop
    for _ in range(num_procs):
        task_queue.put("STOP")

    results.sort(key=lambda tup: tup[0])

    return [item[1] for item in map(list, results)]
