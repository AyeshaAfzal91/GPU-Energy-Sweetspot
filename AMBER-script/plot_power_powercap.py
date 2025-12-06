import pandas as pd
import matplotlib.pyplot as plt
import subprocess

benchmarks = ["Cellulose_NPT", "Cellulose_NVE", "FactorIX_NPT", "FactorIX_NVE", "JAC_NPT", "JAC_NVE", "myoglobin", "nucleosome", "STMV_NPT", "STMV_NVE", "TRPCage"]
architectures = ["A40powercap", "A100powercap"]  # Extend if needed

for arch in architectures:
    for benchmark in benchmarks:
        output_filename = f"{arch}_output_{benchmark}.csv"
        plot_filename = f"{arch}_power_vs_powercap_{benchmark}.png"


        command = f"""
        awk -F', ' '
        BEGIN {{
            OFS = ",";
        }}
        {{
            # Extract powercap from filename: *_<powercap>_powerlog.csv
            match(FILENAME, /_([0-9]+)_powerlog\\.csv$/, arr);
            powercap = arr[1];

            if (powercap == "") next;

            # Extract numeric power value from second field
            split($2, p, " ");
            power = p[1] + 0;

            # Extract utilization from line (regex)
            match($0, /, *([0-9.]+) *%/, u);
            util = u[1] + 0;

            if (util < 85) next;

            power_sum[powercap] += power;
            count[powercap]++;
        }}
        END {{
            print "powercap,avg_power";
            for (pc in power_sum) {{
                avg_power = power_sum[pc] / count[pc];
                print pc "," avg_power;
            }}
        }}
        ' {arch}/*/{benchmark}_*powerlog.csv 2>/dev/null | sort -t',' -k1,1n > {output_filename}
        """
        try:
            subprocess.run(command, shell=True, check=True)
            print(f"Data for {arch} - {benchmark} saved to {output_filename}")
            
            # Load CSV directly
            data = pd.read_csv(output_filename)
            data['powercap'] = pd.to_numeric(data['powercap'], errors='coerce')
            data['avg_power'] = pd.to_numeric(data['avg_power'], errors='coerce')

            print(f"Loaded data shape for {arch} - {benchmark}: {data.shape}")
            print(data.head())

            plt.figure(figsize=(10, 6))

            powercap_vals = data['powercap'].astype(int).to_numpy()
            avg_power_vals = data['avg_power'].to_numpy()

            plt.plot(powercap_vals, avg_power_vals, marker='o', linestyle='-')

            plt.xlabel('Power Cap [W]')
            plt.ylabel('Average Power [W]')
            plt.title(f'Average Power vs. Power Cap for {arch} - {benchmark}')
            plt.xlim(left=0)
            plt.ylim(bottom=0)
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(plot_filename, dpi=300)
            plt.close()
            print(f"Plot for {arch} - {benchmark} saved to {plot_filename}")

        except subprocess.CalledProcessError as e:
            print(f"Error running command for {arch} - {benchmark}: {e}")
        except FileNotFoundError:
            print(f"Error: grep could not find files for {arch} - {benchmark}.")
        except pd.errors.EmptyDataError:
            print(f"Warning: {output_filename} is empty for {arch} - {benchmark}. No plot generated.")

print("Processing and plotting complete for all architectures and benchmarks.")

