"""Module for helping running simulations.
"""
import fcntl


def write_results(results,resultfilename):
    """Safe way of writing results to a file that are used by many processes.
    """
    resultfile=open(resultfilename,"a")
    fcntl.flock(resultfile, fcntl.LOCK_EX)
    for result in results:
        resultfile.write(result)
    fcntl.flock(resultfile, fcntl.LOCK_UN)
    resultfile.close()

