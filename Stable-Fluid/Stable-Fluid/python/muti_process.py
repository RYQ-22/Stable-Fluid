import multiprocessing
import subprocess

def run_script(argument):
    subprocess.run(["python", "script.py", "-input", str(argument)])

def run_scripts_in_parallel(argument, max_process):
    with multiprocessing.Pool(max_process) as pool:
        pool.map(run_script, argument)

arguments = range(100)

run_scripts_in_parallel(arguments, 24)
    