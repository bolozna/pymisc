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


def parse_number(inputstr):
    """Parse numbers given as a string by the user.
    """
    try:
        ilist=[float(inputstr)]
    except ValueError:
        if inputstr.startswith("linspace(") and inputstr[-1]==")":
            inputstr=inputstr.strip("linspace(").strip(")")
            vals=map(float,inputstr.split(","))
            ilist=numpy.linspace(vals[0],vals[1],int(vals[2]))
        elif inputstr.startswith("[") and inputstr[-1]=="]":
            inputstr=inputstr.strip("[").strip("]")
            ilist=map(float,inputstr.split(","))
        else:
            raise Exception("Invalid parameter: "+inputstr)
    return ilist
