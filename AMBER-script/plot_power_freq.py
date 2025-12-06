import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os

benchmarks = ["Cellulose_NPT", "Cellulose_NVE", "FactorIX_NPT", "FactorIX_NVE", "JAC_NPT", "JAC_NVE", "myoglobin", "nucleosome", "STMV_NPT", "STMV_NVE", "TRPCage"]
#architectures = ["A40PowerCap", "A100PowerCap"]  # Extend if needed
architectures = ["A100freq", "A40freq"]
for arch in architectures:
    for benchmark in benchmarks:
        output_filename = f"{arch}_output_{benchmark}.csv"
        plot_filename = f"{arch}_power_vs_freq2_{benchmark}.png"


        command = f"""
        awk '
        /W/ {{
            # Extract frequencies from filename (comma OR dash)
            match(FILENAME, /_([0-9]+)[,-]([0-9]+)_powerlog\\.csv$/, arr);
            freq1 = arr[1];
            freq2 = arr[2];

            if (freq1 == "" || freq2 == "") next;

            # Extract power from line (regex)
            match($0, /, *([0-9.]+) W/, p);
            power = p[1] + 0;

            # Extract utilization (regex)
            match($0, /, *([0-9.]+) *%/, u);
            util = u[1] + 0;

            if (util <= 85) next;

            key = freq1 "," freq2;
            power_sum[key] += power;
            count[key]++;
        }}
        END {{
            print "freq1,freq2,avg_power";
            for (k in power_sum) {{
                avg_power = power_sum[k] / count[k];
                split(k, a, ",");
                print a[1] "," a[2] "," avg_power;
            }}
        }}
            ' {arch}/*/{benchmark}_*powerlog.csv 2>/dev/null | sort -t',' -k2,2n > {output_filename}
        """

        try:
            subprocess.run(command, shell=True, check=True)
            print(f"Data for {arch} - {benchmark} saved to {output_filename}")

            try:
                data = pd.read_csv(output_filename)
                #data.columns = data.iloc[0]   # first row becomes header
                #data = data[1:]               # drop header row
                #data = data.reset_index(drop=True)
                print(f"Loaded data shape for {arch} - {benchmark}: {data.shape}")
                print(data.head())

                plt.figure(figsize=(10, 6))

                for freq1_val in sorted(data['freq1'].unique()):
                    subset = data[data['freq1'] == freq1_val]
                    freq2_np = subset['freq2'].to_numpy()
                    avg_power_np = subset['avg_power'].to_numpy()
                    plt.plot(freq2_np, avg_power_np,
                             marker='o', linestyle='-', label=f'Memory frequency [MHz]={freq1_val}')

                plt.xlabel('Graphics frequency [MHz]')
                plt.ylabel('Average Power [W]')
                plt.title(f'Average Power vs. Graphics frequency for {arch} - {benchmark}')
                plt.xlim(left=0)
                plt.ylim(bottom=0)
                plt.grid(True)
                plt.legend()
                plt.tight_layout()
                plt.savefig(plot_filename, dpi=300)
                plt.close()
                print(f"Plot for {arch} - {benchmark} saved to {plot_filename}")

            except FileNotFoundError:
                print(f"Error: {output_filename} not found for plotting.")
            except pd.errors.EmptyDataError:
                print(f"Warning: {output_filename} is empty for {arch} - {benchmark}. No plot generated.")

        except subprocess.CalledProcessError as e:
            print(f"Error running command for {arch} - {benchmark}: {e}")
        except FileNotFoundError:
            print(f"Error: grep could not find files for {arch} - {benchmark}.")

print("Processing and plotting complete for all architectures and benchmarks.")

