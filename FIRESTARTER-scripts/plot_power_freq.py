#!/usr/bin/env python3

import csv
import matplotlib.pyplot as plt
import os
import re
import sys

class Result:
    def __init__(self, powerdraw, duration, gflops):
        self.powerdraw = powerdraw
        self.duration = duration
        self.gflops = gflops

def scan_results():
    # Iterate over all directories and search those, which match pattern "FIRESTARTER-*-{frequency,powercap}"
    for d in os.scandir(results_dir):
        if not d.is_dir():
            continue

        f_match = re.fullmatch("^FIRESTARTER-(\\w+)-frequency$", os.path.basename(d.path))
        if f_match:
            arch = f_match[1]
            if arch not in frequency_results:
                frequency_results[arch] = {}
            scan_subresults(d.path, arch, True)
            continue
        p_match = re.fullmatch("^FIRESTARTER-(\\w+)-powercap$", os.path.basename(d.path))
        if p_match:
            arch = p_match[1]
            if arch not in powercap_results:
                powercap_results[arch] = {}
            scan_subresults(d.path, arch, False)
            continue
        # You can add more variants here

def scan_subresults(path, arch, is_freq):
    # Iterate over all subdirectories, which look like a number
    for d in os.scandir(path):
        if not d.is_dir():
            continue
        d_match = re.fullmatch("^(\\d+)$", os.path.basename(d.path))
        if not d_match:
            continue

        csv_path = scan_nvidiasmi_csv(d.path)
        log_path = scan_firestarter_log(d.path)

        parse_run(csv_path, log_path, arch, is_freq, int(d_match[1]))

def scan_nvidiasmi_csv(path, arch, is_freq, value):
    # Iterate over all files, which looks like *.csv
    # Take the first one and parse that one
    for f in os.scandir(path):
        if not f.is_file():
            continue
        f_match = re.fullmatch("^.*\\.csv$", os.path.basename(f.path))
        if not f_match:
            continue
        return f.path
    raise Exception("Unable to find csv file in directory: {}".format(path))

def scan_firestarter_log(path):
    # Iterate over all files, which looks like *.fsout
    # Take the first one and parse that one
    for f in os.scandir(path):
        if not f.is_file():
            continue
        f_match = re.fullmatch("^.*\\.fsout$", os.path.basename(f.path))
        if not f_match:
            continue
        return f.path
    raise Exception("Unable to find csv file in directory: {}".format(path))

def parse_run(csv_path, log_path, arch, is_freq, value):
    # Parse nvidia-smi csv
    power_sum = 0.0
    power_iters = 0
    warmup_iters = 20
    with open(csv_path, "r") as f:
        rd = csv.reader(f)
        _ = next(rd) # skip header

        for row in rd:
            # parse "9 %" to 9
            gpu_util = float(row[2].strip().split(" ")[0])
            # parse "12.34 W" to 12.34
            power_draw = float(row[1].strip().split(" ")[0])
            if gpu_util < 98.0:
                continue

            # ignore the first 20 results, since the power readings will be off due to moving average filter
            if warmup_iters > 0:
                warmup_iters -= 1
                continue

            power_sum += power_draw
            power_iters += 1

    if power_iters == 0:
        raise Exception("Not enough results in csv: {}", csv_path)

    # TODO fix this
    if is_freq:
        frequency_results[arch][value] = power_sum / power_iters
    else:
        powercap_results[arch][value] = power_sum / power_iters

    # TODO parse GFLOPS from firestarter log:
    with open(log_path, "r") as log_file:
        for line in log_file:
            if 
    # ^.*: ((?:0|[1-9]\d*)(?:\.\d+)?) GFLOPS.*$

def create_power_plots():
    results = [
        {
            "name": "powercap",
            "xlabel": "Power cap (W)",
            "vs": "Power cap",
            "result": powercap_results,
        },
        {
            "name": "frequency",
            "xlabel": "Graphics frequency (MHz)",
            "vs": "Graphics frequency",
            "result": frequency_results,
        },
    ]

    for result in results:
        plt.figure(figsize=(9, 6))
        plt.xlabel(result["xlabel"])
        plt.ylabel("Average GPU power draw (W)")
        plt.title("Avg GPU power draw vs {}".format(result["vs"]))
        plt.legend(fontsize=8, ncol=2)
        plt.grid(True)

        for arch_name, arch_dict in result["result"].items():
            keys = sorted(arch_dict)
            values = [arch_dict[k] for k in keys]
            plt.plot(keys, values, marker="o", linestyle="-", linewidth=1, markersize=3, label=arch_name)

        plt.ylim(ymin=0)
        plt.xlim(xmin=0)
        plt.legend()

        figure_dir = os.path.join(results_dir, "FIRESTARTER-plots")
        figure_path = os.path.join(figure_dir, "{}.svg".format(result["name"]))
        os.makedirs(figure_dir, exist_ok=True)
        plt.savefig(figure_path, dpi=300, bbox_inches="tight")
        plt.close()

def create_efficiency_plots():
    for result in results:
        pass

# Get the parent directory of this script
results_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

# *_results["A40"][210] = 234
powercap_results = {}
frequency_results = {}

scan_results()
create_power_plots()
