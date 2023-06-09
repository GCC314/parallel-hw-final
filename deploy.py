import subprocess

PERF_DATA_DIR = "/root/perf"

subprocess.run("./utils/init.sh",shell=True)
subprocess.run("make clean",shell=True)
subprocess.run("make pre",shell=True)
subprocess.run("make lapack",shell=True)
subprocess.run("make scalapack",shell=True)
subprocess.run("make main",shell=True)

def runperf(args):
    v = args[0]
    points = args[1]
    thread_cnt = args[2]
    test_rootdir = PERF_DATA_DIR
    test_dir = test_rootdir
    input_file = "".join([
        f"isHexahedral 0\n",
        f"lx 1000\n",
        f"ly 1000\n",
        f"lz 1000\n",
        f"thetaxy 0\n",
        f"thetayz 0\n",
        f"thetaxz 0\n",
        f"support_SH 0\n",
        f"diago_lib scalapack\n",
        f"support_Periodic_Boundary 0\n",
        f"multi_parallel_strategies 0\n",
        f"points_path {test_dir}/points{points}.txt\n",
        f"v_path {test_dir}/V{v}.txt\n",
        f"distribution_path {test_dir}/distribution80.txt"
    ])
    with open("INPUT.txt", "w") as fp: print(input_file, file=fp)
    subprocess.run(f"make run_with_memcheck OPENMP_THREAD_NUM={thread_cnt} INPUT_FILE_DIR=INPUT.txt &> output/{v}-{points}-{thread_cnt}.log",shell=True)
    print(args, "done")


for V in [128, 256, 512, 1024]:
    for p in [2, 10, 50]:
        for th in [1, 8, 16, 32]:
            runperf((V, p, th))