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
Parallelizer.py

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
    """
    Abstract parallelization class
    """

    def __init__(self, mode=None, num_procs=None):
        """
        This will initialize the Parallelizer class and kick off the specific classes for multiprocessing.

        Default num_procs is all the processesors possible

        This will also establish:
            :self   str     self.mode:  The mode which will be used for this paralellization. This determined by mode and the enviorment
                                        default is multiprocessing on all processors unless stated otherwise.
            :self   class   self.parallel_obj: This is the obstantiated object of the class of parallizations.
                            ie)
                            self.parallel_obj will be set to None for simpler parallization methods like serial
            :self   int     self.num_processor:   the number of processors or nodes that will be used. If None than we will use all available nodes/processors
                                                This will be overriden and fixed to a single processor if mode==serial
        Inputs:
        :param str mode: the multiprocess mode to be used, ie) serial, multiprocessing, or None:
                            if None then we will try to pick a possible multiprocessing choice. This should only be used for
                            top level coding. It is best practice to specify which multiprocessing choice to use.
        :param int num_procs:   the number of processors or nodes that will be used. If None than we will use all available nodes/processors
                                        This will be overriden and fixed to a single processor if mode==serial
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
            self.num_procs = self.compute_nodes()

        else:
            self.num_procs = num_procs

    def run(self, args, func, num_procs=None, mode=None):
        """
        Run a task in parallel across the system.

        Mode can be one of 'multiprocessing' or 'none' (serial). If it is not
        set, the best value will be determined automatically.

        By default, this method will use the full resources of the system. However,
        if the mode is set to 'multiprocessing', num_procs can control the number
        of threads initialized when it is set to a nonzero value.

        Example: If one wants to multiprocess function  def foo(x,y) which takes 2 ints and one wants to test all permutations of x and y between 0 and 2:
                    args = [(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)]
                    func = foo      The namespace of foo


        Inputs:
        :param python_obj func: This is the object of the function which will be used.
        :param list args: a list of lists/tuples, each sublist/tuple must contain all information required by the function for a single object which will be multiprocessed
        :param int num_procs:  (Primarily for Developers)  the number of processors or nodes that will be used. If None than we will use all available nodes/processors
                                        This will be overriden and fixed to a single processor if mode==serial
        :param str mode:  (Primarily for Developers) the multiprocess mode to be used, ie) serial, multiprocessing, or None:
                            if None then we will try to pick a possible multiprocessing choice. This should only be used for
                            top level coding. It is best practice to specify which multiprocessing choice to use.
                            if you have smaller programs used by a larger program, with both mpi enabled there will be problems, so specify multiprocessing is important.
                            BEST TO LEAVE THIS BLANK. NOTE: mpi no longer supported anyway.
        Returns:
        :returns: list results: A list containing all the results from the multiprocess
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
            return MultiThreading(args, num_procs, func)
        else:
            # serial is running the ParallelThreading with num_procs=1
            return MultiThreading(args, 1, func)

    def compute_nodes(self, mode=None):
        """
        Computes the number of "compute nodes" according to the selected mode.

        For multiprocessing this is the number of available cores
        For serial, this value is 1
        Returns:
        :returns: int num_procs: the number of nodes/processors which is to be used
        """
        if mode is None:
            mode = self.mode

        return multiprocessing.cpu_count() if mode == "multiprocessing" else 1

    def return_node(self):
        """
        Returns the number of "compute nodes" according to the selected mode.

        For multiprocessing this is the number of available cores
        For serial, this value is 1
        Returns:
        :returns: int num_procs: the number of nodes/processors which is to be used
        """
        return self.num_procs


def MultiThreading(inputs, num_procs, task_name):
    """Initialize this object.

    Args:
        inputs ([data]): A list of data. Each datum contains the details to
            run a single job on a single processor.
        num_procs (int): The number of processors to use.
        task_class_name (class): The class that governs what to do for each
            job on each processor.
    """

    results = []

    # If there are no inputs, just return an empty list.
    if len(inputs) == 0:
        return results

    inputs = check_and_format_inputs_to_list_of_tuples(inputs)

    num_procs = count_processors(len(inputs), num_procs)

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
        results = start_processes(tasks, num_procs)

    return results


###
# Worker function
###


def worker(input, output):
    for seq, job in iter(input.get, "STOP"):
        func, args = job
        result = func(*args)
        ret_val = (seq, result)
        output.put(ret_val)


def check_and_format_inputs_to_list_of_tuples(args):
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


def count_processors(num_inputs, num_procs):
    """
    Checks processors available and returns a safe number of them to
    utilize.

    :param int num_inputs: The number of inputs.
    :param int num_procs: The number of desired processors.

    :returns: The number of processors to use.
    """
    # first, if num_procs <= 0, determine the number of processors to
    # use programatically
    if num_procs <= 0:
        num_procs = multiprocessing.cpu_count()

    # reduce the number of processors if too many have been specified
    if num_inputs < num_procs:
        num_procs = num_inputs

    return num_procs


def start_processes(inputs, num_procs):
    """
    Creates a queue of inputs and outputs
    """

    # Create queues
    task_queue = multiprocessing.Queue()
    done_queue = multiprocessing.Queue()

    # Submit tasks
    for item in inputs:
        task_queue.put(item)

    # Start worker processes
    for _ in range(num_procs):
        multiprocessing.Process(target=worker, args=(task_queue, done_queue)).start()

    results = [done_queue.get() for _ in range(len(inputs))]
    # Tell child processes to stop
    for _ in range(num_procs):
        task_queue.put("STOP")

    results.sort(key=lambda tup: tup[0])

    return [item[1] for item in map(list, results)]
