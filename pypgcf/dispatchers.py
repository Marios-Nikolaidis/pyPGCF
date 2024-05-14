from os import fork, system
from subprocess import Popen
from concurrent.futures import ProcessPoolExecutor
from typing import Callable, List, Union
from tqdm import tqdm

def multiprocess_dispatch(f: Union[Callable,str], args: list, num_procs: int, show_progress: bool, description: str = "") -> List:
    """
    Map a list of executable to a ProcessPoolExecutor
    Return: list
    """
    if isinstance(f, str): # If it is a string then just call the system
        f = system
    with ProcessPoolExecutor(num_procs) as executor:
        if show_progress:
            results = list(
                tqdm(
                    executor.map(f, args),
                    total=len(args),
                    desc=description,
                    ascii=True,
                    leave=True,
                )
            )
        else:
            results = list(executor.map(f, args))
    return results

def execute_command(cmd: str) -> int:
    ret = system(cmd)
    return ret

