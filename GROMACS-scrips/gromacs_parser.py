import os
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.optimize import curve_fit

# ==========================================
# CONFIGURATION
# ==========================================
ROOT = "."
ARCHS = ["gromacs_H100_cap", "gromacs_H200_cap", "gromacs_H100_freq", "gromacs_H200_freq"]
UTIL_THRESHOLD = 80.0
OUTDIR = "plots_gpu"
CSV_CACHE = "gpu_gromacs_data.csv"

ATOM_COUNTS = {
    "FactorIX_NVE": 90906, "FactorIX_NPT": 90906,
    "JAC_NVE": 23558, "JAC_NPT": 23558,
    "Cellulose_NVE": 408609, "Cellulose_NPT": 408609,
    "STMV_NVE": 1067095, "STMV_NPT": 1067095,
    "nucleosome": 25095, "TRPCage": 304, "myoglobin": 2492,
    "2md_start0": 20248, "FL_md1_berendsen": 170320, 
    "rnanvt": 31889, "eag1": 615924, "PI_large_test": 80289, "stmv_pme_nvt": 1066628
}

# --- UNIFIED REGEXES ---
FREQ_RE = re.compile(r"^(?P<bench>.+?)_(?P<mem>\d+)-(?P<gfx>\d+)_.*powerlog\.csv$", re.IGNORECASE)
CAP_RE = re.compile(r"^(?P<bench>.+?)_(?P<pc>\d+)_.*powerlog\.csv$", re.IGNORECASE)

# ==========================================
# HELPER FUNCTIONS
# ==========================================
def avg_power_from_powerlog(path, util_thr=85.0):
    vals = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.lower().startswith("timestamp"): continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 3: continue
            mp, mu = re.search(r"([0-9.]+)", parts[1]), re.search(r"([0-9.]+)", parts[2])
            if not mp or not mu: continue
            p, u = float(mp.group(1)), float(mu.group(1))
            if u >= util_thr: vals.append(p)
    return (sum(vals) / len(vals)) if vals else None

def find_gromacs_log(run_dir, bench, search_suffix):
    search_str = f"{bench}_{search_suffix}".lower()
    for fn in os.listdir(run_dir):
        low = fn.lower()
        if (low.endswith(".out") or low.endswith(".log")) and search_str in low:
            return os.path.join(run_dir, fn)
    return None

def parse_nsday_from_gromacs(log_path):
    with open(log_path, "r", errors="ignore") as f: txt = f.read()
    m = re.search(r"Performance:\s+([0-9.]+)", txt)
    if m: return float(m.group(1))
    return None

def iqr_filter(df, col, k=1.5):
    q1, q3 = df[col].quantile(0.25), df[col].quantile(0.75)
    iqr = q3 - q1
    return df[(df[col] >= q1 - k * iqr) & (df[col] <= q3 + k * iqr)]

def get_atom_count(bench_name):
    if bench_name in ATOM_COUNTS: return ATOM_COUNTS[bench_name]
    for key, val in ATOM_COUNTS.items():
        if bench_name.lower() == key.lower(): return val
    return None

# ==========================================
# PARSING & CACHING LOGIC
# ==========================================
def load_or_parse_data(clean_cache):
    if clean_cache and os.path.exists(CSV_CACHE):
        print(f"🧹 [--clean] Flag detected. Deleting cached {CSV_CACHE}...")
        os.remove(CSV_CACHE)

    if os.path.exists(CSV_CACHE):
        print(f"📦 Loading cached data from {CSV_CACHE} (Use --clean to re-parse logs)")
        return pd.read_csv(CSV_CACHE)

    print("🔍 Parsing raw GROMACS logs. This might take a moment...")
    rows = []
    for arch in ARCHS:
        arch_dir = os.path.join(ROOT, arch)
        if not os.path.isdir(arch_dir): continue

        for run_root, _, files in os.walk(arch_dir):
            for fn in files:
                run_type, bench, search_suffix = None, None, None
                mem_mhz, gfx_mhz, pc_w = pd.NA, pd.NA, pd.NA
                
                m_freq = FREQ_RE.match(fn)
                if m_freq:
                    run_type, bench = "frequency", m_freq.group("bench")
                    mem_mhz, gfx_mhz = int(m_freq.group("mem")), int(m_freq.group("gfx"))
                    search_suffix = f"{mem_mhz}-{gfx_mhz}"
                else:
                    m_cap = CAP_RE.match(fn)
                    if m_cap:
                        run_type, bench = "powercap", m_cap.group("bench")
                        pc_w = int(m_cap.group("pc"))
                        search_suffix = f"{pc_w}"
                    else:
                        continue

                powerlog = os.path.join(run_root, fn)
                avg_p = avg_power_from_powerlog(powerlog, UTIL_THRESHOLD)
                if avg_p is None: continue

                gmx_log = find_gromacs_log(run_root, bench, search_suffix)
                if gmx_log is None: continue

                ns_day = parse_nsday_from_gromacs(gmx_log)
                if ns_day is None: continue

                rows.append(dict(
                    arch=arch, run_type=run_type, benchmark=bench,
                    mem_freq_mhz=mem_mhz, gfx_freq_mhz=gfx_mhz, powercap_w=pc_w,
                    ns_day=ns_day, avg_power_w=avg_p, run_dir=run_root, 
                    log_file=gmx_log, powerlog=powerlog
                ))

    df = pd.DataFrame(rows)
    if df.empty: raise SystemExit("No data parsed. Check log files and directories.")

    df["atom_count"] = df["benchmark"].apply(get_atom_count)
    df = df.dropna(subset=["atom_count"]).copy()
    df["atom_count"] = df["atom_count"].astype(int)

    # Calculate metrics before caching
    df["atom_ns_day"] = df["ns_day"] * df["atom_count"]
    df["efficiency"] = df["atom_ns_day"] / df["avg_power_w"]
    df["edp"] = df["avg_power_w"] / (df["atom_ns_day"] ** 2)

    df.to_csv(CSV_CACHE, index=False)
    print(f"✅ Data parsed and cached to {CSV_CACHE}")
    return df

# ==========================================
# MAIN SCRIPT (CLI & MODULAR PLOTTING)
# ==========================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modular GPU Benchmark Plotter for GROMACS")
    
    parser.add_argument("--clean", action="store_true", help="Force re-parsing of raw logs")
    parser.add_argument("--all", action="store_true", help="Generate ALL plots and reports")
    parser.add_argument("--tradeoff", action="store_true", help="Generate Trade-off Analysis CSV")
    parser.add_argument("--z-normal", action="store_true", help="Generate Normal Z-plots (Size-weighted Performance vs Power)")
    parser.add_argument("--z-energy", action="store_true", help="Generate Publication Energy Z-plots")
    parser.add_argument("--power-freq", action="store_true", help="Generate Avg Power vs Graphics Frequency plots (Frequency runs only)")
    parser.add_argument("--power-cap", action="store_true", help="Generate Avg Power vs Power Cap plots (Powercap runs only)")
    parser.add_argument("--model-power", action="store_true", help="Generate Predictive Power Models (Levenberg-Marquardt) for Frequency runs")

    args = parser.parse_args()

    if not any([args.all, args.tradeoff, args.z_normal, args.z_energy, args.power_freq, args.power_cap, args.model_power]):
        print("⚠️ No output flags specified. Only parsing/caching data. Use --help to see available options.")

    os.makedirs(OUTDIR, exist_ok=True)
    df = load_or_parse_data(args.clean)

    # Split datasets early for the specific plots
    df_freq = df[df["run_type"] == "frequency"]
    df_cap = df[df["run_type"] == "powercap"]

    # ----------------------------
    # 0) Trade-off Analysis
    # ----------------------------
    if args.all or args.tradeoff:
        print("Generating Trade-off Analysis...")
        tradeoff_rows = []
        for arch, sub_arch in df.groupby("arch"):
            for bench, ssub in sub_arch.groupby("benchmark"):
                if ssub.empty: continue
                
                # 1. Identify the three key operating points
                max_perf_row = ssub.loc[ssub["atom_ns_day"].idxmax()]
                opt_eff_row = ssub.loc[ssub["efficiency"].idxmax()]
                min_edp_row = ssub.loc[ssub["edp"].idxmin()]

                # 2. Extract Performance and Efficiency for each point
                max_perf = max_perf_row["atom_ns_day"]
                max_perf_eff = max_perf_row["efficiency"]

                eff_perf = opt_eff_row["atom_ns_day"]
                eff_eff = opt_eff_row["efficiency"]
                
                edp_perf = min_edp_row["atom_ns_day"]
                edp_eff = min_edp_row["efficiency"]

                # 3. Calculate Trade-offs relative to Max Performance
                
                # (A) Trade-off when targeting Maximum Efficiency
                perf_drop_to_eff_pct = ((max_perf - eff_perf) / max_perf) * 100
                eff_gain_to_eff_pct = ((eff_eff - max_perf_eff) / max_perf_eff) * 100

                # (B) Trade-off when targeting Minimum EDP (Sweet Spot)
                perf_drop_to_edp_pct = ((max_perf - edp_perf) / max_perf) * 100
                eff_gain_to_edp_pct = ((edp_eff - max_perf_eff) / max_perf_eff) * 100

                tradeoff_rows.append({
                    "Architecture": arch, 
                    "Run_Type": ssub["run_type"].iloc[0], 
                    "Benchmark": bench,
                    "Max_Perf (atom-ns/day)": max_perf, 
                    "Max_Perf_Efficiency": max_perf_eff,
                    
                    "Opt_Eff_Perf (atom-ns/day)": eff_perf, 
                    "Opt_Eff_Efficiency": eff_eff,
                    "Perf_Drop_to_Max_Eff (%)": perf_drop_to_eff_pct, 
                    "Efficiency_Gain_to_Max_Eff (%)": eff_gain_to_eff_pct,
                    
                    "Min_EDP_Perf (atom-ns/day)": edp_perf,
                    "Min_EDP_Efficiency": edp_eff,
                    "Perf_Drop_to_Min_EDP (%)": perf_drop_to_edp_pct,
                    "Efficiency_Gain_to_Min_EDP (%)": eff_gain_to_edp_pct
                })

        tradeoff_df = pd.DataFrame(tradeoff_rows)
        out_csv = os.path.join(OUTDIR, "efficiency_tradeoff_gromacs.csv")
        tradeoff_df.to_csv(out_csv, index=False)
        print(f" -> Saved to {out_csv}")

    # ----------------------------
    # 1) Normal Z-plot (Performance vs Power)
    # ----------------------------
    if args.all or args.z_normal:
        print("Generating Normal Z-plots...")
        for arch, sub in df.groupby("arch"):
            plt.figure(figsize=(9, 6))
            for bench, ssub in sub.groupby("benchmark"):
                if ssub["run_type"].iloc[0] == "frequency":
                    ssub = ssub.sort_values(by="gfx_freq_mhz")
                else:
                    ssub = ssub.sort_values(by="powercap_w")

                ssub = iqr_filter(ssub, "atom_ns_day", k=1.0)
                if ssub.empty: continue

                plt.scatter(ssub["atom_ns_day"], ssub["avg_power_w"], s=12, label=bench)
                plt.plot(ssub["atom_ns_day"], ssub["avg_power_w"], linewidth=0.8)

            plt.xlabel("Size-weighted performance (atom-ns/day)")
            plt.ylabel(f"GPU Avg Power (W) (util ≥ {UTIL_THRESHOLD:.0f}%)")
            plt.title(f"GPU Z-plot: Power vs Size-weighted Performance — {arch}")
            plt.legend(fontsize=8, ncol=2)
            plt.savefig(os.path.join(OUTDIR, f"{arch}_gpu_zplot_atom_weighted.png"), dpi=300, bbox_inches="tight")
            plt.close()

    # ----------------------------
    # 2) Publication Energy Z-plot
    # ----------------------------
    if args.all or args.z_energy:
        print("Generating Publication Energy Z-plots...")
        for arch, sub in df.groupby("arch"):
            plt.figure(figsize=(9, 6))
            for bench, ssub in sub.groupby("benchmark"):
                is_freq = ssub["run_type"].iloc[0] == "frequency"
                if is_freq:
                    ssub = ssub.sort_values(by="gfx_freq_mhz")
                else:
                    ssub = ssub.sort_values(by="powercap_w")

                ssub = iqr_filter(ssub, "atom_ns_day", k=1.0)
                if ssub.empty: continue

                plt.scatter(ssub["atom_ns_day"], ssub["efficiency"], s=12, label=bench)
                line_color = plt.plot(ssub["atom_ns_day"], ssub["efficiency"], linewidth=0.8)[0].get_color()

                max_eff_idx, sweet_spot_idx = ssub["efficiency"].idxmax(), ssub["edp"].idxmin()
                eff_row, sweet_row = ssub.loc[max_eff_idx], ssub.loc[sweet_spot_idx]

                eff_perf, eff_eff = eff_row["atom_ns_day"], eff_row["efficiency"]
                eff_val = int(eff_row['gfx_freq_mhz']) if is_freq else int(eff_row['powercap_w'])
                eff_label_text = f"{eff_val} {'MHz' if is_freq else 'W'}"

                opt_perf, opt_eff = sweet_row["atom_ns_day"], sweet_row["efficiency"]
                opt_val = int(sweet_row['gfx_freq_mhz']) if is_freq else int(sweet_row['powercap_w'])
                edp_label_text = f"{opt_val} {'MHz' if is_freq else 'W'}"

                if max_eff_idx == sweet_spot_idx:
                    plt.scatter(eff_perf, eff_eff, color='black', s=45, marker="D", zorder=5)
                    plt.annotate(eff_label_text, (eff_perf, eff_eff), textcoords="offset points", 
                                 xytext=(0, 12), ha='center', fontsize=8, rotation=0, color=line_color)
                    #plt.annotate(edp_label_text, (opt_perf, opt_eff), textcoords="offset points", 
                                 #xytext=(10, 10), ha='left', fontsize=8, rotation=-45, color=line_color)
                else:
                    plt.scatter(eff_perf, eff_eff, color='black', s=40, marker="s", zorder=5)
                    plt.annotate(eff_label_text, (eff_perf, eff_eff), textcoords="offset points", 
                                 xytext=(0, 8), ha='center', fontsize=8, rotation=0, color=line_color)

                    plt.scatter(opt_perf, opt_eff, color='black', s=40, marker="o", zorder=5)
                    plt.annotate(edp_label_text, (opt_perf, opt_eff), textcoords="offset points", 
                                 xytext=(-8, -8), ha='left', fontsize=8, rotation=-45, color=line_color)

            plt.xlabel("Size-weighted performance (atom-ns/day)")
            plt.ylabel(f"Efficiency (atom-ns/day/W)")
            #plt.title(f"Energy Z-plot — {arch}")
            
            # Custom Mixed Legend
            handles, labels = plt.gca().get_legend_handles_labels()
            marker_max = Line2D([0], [0], marker='s', color='w', markerfacecolor='black', markersize=7, label='Max $\\eta$')
            marker_edp = Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=7, label='Min EDP')
            marker_both = Line2D([0], [0], marker='D', color='w', markerfacecolor='black', markersize=7, label='Coincide')
            handles.extend([marker_max, marker_edp, marker_both])
            
            plt.legend(handles=handles, fontsize=8, ncol=2, loc='lower right')
            plt.grid(True, linestyle="--", alpha=0.4)
            plt.savefig(os.path.join(OUTDIR, f"{arch}_gpu_zplot_energy.png"), dpi=300, bbox_inches="tight")
            plt.close()

    # ----------------------------
    # 3) Power vs Frequency (df_freq only)
    # ----------------------------
    if args.all or args.power_freq:
        print("Generating Power vs Frequency plots...")
        if not df_freq.empty:
            for arch, sub_arch in df_freq.groupby("arch"):
                plt.figure(figsize=(9, 6))
                for bench, ssub in sub_arch.groupby("benchmark"):
                    mem_choice = int(ssub["mem_freq_mhz"].mode().iloc[0])
                    s = ssub[ssub["mem_freq_mhz"] == mem_choice].sort_values(by="gfx_freq_mhz")
                    if s.empty: continue
                    plt.plot(s["gfx_freq_mhz"], s["avg_power_w"], marker="o", linestyle="-", linewidth=1, markersize=3, label=bench)
                
                plt.xlabel("Graphics frequency (MHz)")
                plt.ylabel("Average GPU Power (W)")
                plt.title(f"Avg GPU Power vs Graphics frequency — {arch}")
                plt.legend(fontsize=8, ncol=2)
                plt.grid(True)
                plt.ylim(ymin=0)
                plt.savefig(os.path.join(OUTDIR, f"{arch}_avgpower_vs_gfxfreq.png"), dpi=300, bbox_inches="tight")
                plt.close()
        else:
            print(" -> Skipping: No frequency data found.")

    # ----------------------------
    # 4) Power vs Powercap (df_cap only)
    # ----------------------------
    if args.all or args.power_cap:
        print("Generating Power vs Power Cap plots...")
        if not df_cap.empty:
            for arch, sub_arch in df_cap.groupby("arch"):
                plt.figure(figsize=(9, 6))
                for bench, ssub in sub_arch.groupby("benchmark"):
                    s = ssub.sort_values(by="powercap_w")
                    if s.empty: continue
                    plt.plot(s["powercap_w"], s["avg_power_w"], marker="o", linestyle="-", linewidth=1, markersize=3, label=bench)

                plt.xlabel("Power cap (W)")
                plt.ylabel("Average GPU Power (W)")
                plt.title(f"Avg GPU Power vs Power cap — {arch}")
                plt.legend(fontsize=8, ncol=2)
                plt.grid(True)
                plt.savefig(os.path.join(OUTDIR, f"{arch}_avgpower_vs_powercap.png"), dpi=300, bbox_inches="tight")
                plt.close()
        else:
            print(" -> Skipping: No powercap data found.")


    # ----------------------------
    # 6) Piecewise Power Modeling (Linear -> Capped Exponential)
    # ----------------------------
    if args.all or args.model_power:
        print("Generating Piecewise Power Models (Linear -> Capped Exponential)...")
        
        df_freq = df[df["run_type"] == "frequency"]
        
        if not df_freq.empty:
            for arch, sub_arch in df_freq.groupby("arch"):
                plt.figure(figsize=(11, 7))
                custom_handles = []

                for bench, ssub in sub_arch.groupby("benchmark"):
                    ssub = ssub.dropna(subset=["gfx_freq_mhz", "avg_power_w"]).sort_values(by="gfx_freq_mhz")
                    
                    if len(ssub) < 6:
                        print(f" -> Skipping model for {bench}: Need at least 6 points for piecewise fitting.")
                        continue
                    
                    x_data = ssub["gfx_freq_mhz"].to_numpy()
                    y_data = ssub["avg_power_w"].to_numpy()

                    # Plot the raw empirical data points
                    p = plt.plot(x_data, y_data, marker="o", linestyle="", alpha=0.5)
                    color = p[0].get_color()

                    # Find the physical hardware saturation limit (max recorded power)
                    # We add a tiny 0.1W buffer so the math doesn't clip the highest empirical dot awkwardly
                    pwr_limit = y_data.max() + 0.1

                    # Define the two mathematical models
                    def lin_model(x, m, c):
                        return m * x + c
                        
                    def poly_capped_model(x, a, b, c):
                        # Calculate the raw quadratic polynomial growth
                        raw_poly = a * (x**2) + b * x + c
                        # Cap the output at the hardware's maximum power limit
                        return raw_poly

                    best_sse = float('inf')
                    best_split_idx = -1
                    best_lin_popt, best_poly_popt = None, None
                    
                    # Brute force search for the best inflection point
                    for i in range(3, len(x_data) - 2):
                        x_lin, y_lin = x_data[:i], y_data[:i]
                        x_poly, y_poly = x_data[i-1:], y_data[i-1:] 
                        
                        try:
                            lin_popt, _ = curve_fit(lin_model, x_lin, y_lin)
                            sse_lin = np.sum((y_lin - lin_model(x_lin, *lin_popt))**2)
                            
                            # Polynomial initial guesses: a (tiny), b (small), c (baseline power)
                            poly_popt, _ = curve_fit(poly_capped_model, x_poly, y_poly, p0=[1e-5, 1e-2, min(y_poly)], maxfev=5000)
                            sse_poly = np.sum((y_poly - poly_capped_model(x_poly, *poly_popt))**2)
                            
                            if sse_lin + sse_poly < best_sse:
                                best_sse = sse_lin + sse_poly
                                best_split_idx = i
                                best_lin_popt = lin_popt
                                best_poly_popt = poly_popt
                        except RuntimeError:
                            continue

                    if best_split_idx != -1:
                        split_freq = x_data[best_split_idx - 1]
                        
                        # Extract the optimized mathematical parameters
                        m, intercept = best_lin_popt
                        a, b, c = best_poly_popt
                        
                        x_smooth_lin = np.linspace(x_data.min(), split_freq, 50)
                        x_smooth_poly = np.linspace(split_freq, x_data.max(), 100)
                        
                        # Plot the solid linear line and dashed polynomial line
                        plt.plot(x_smooth_lin, lin_model(x_smooth_lin, *best_lin_popt), linestyle="--", color=color, linewidth=2)
                        plt.plot(x_smooth_poly, poly_capped_model(x_smooth_poly, *best_poly_popt), linestyle="--", color=color, linewidth=2)
                        
                        # Calculate Deviation (RMSE)
                        rmse = np.sqrt(best_sse / len(x_data))
                        
                        # Format the equations to clearly show the min() limit
                        # We use scientific notation (.2e) for a and b because they are usually very small decimals
                        lin_eq = f"Lin (f \u2264 {split_freq}): P = {m:.3f}f {intercept:+.1f}"
                        poly_eq = f"Poly (f > {split_freq}): P = min({a:.2e}f\u00b2 {b:+.2e}f {c:+.1f}, {pwr_limit:.1f})"
                        err_text = f"Deviation: \u00B1{rmse:.2f} W"
                        
                        label_text = f"{bench}\n  {lin_eq}\n  {poly_eq}\n  {err_text}"
                        custom_handles.append(Line2D([0], [0], color=color, linestyle="-", marker="o", label=label_text))
                    else:
                        print(f" -> Warning: Piecewise fit failed to converge for {bench}")

                plt.xlabel("Graphics frequency (MHz)")
                plt.ylabel("Average GPU Power (W)")
                plt.title(f"Piecewise Power Modeling — {arch}")
                
                if custom_handles:
                    plt.legend(handles=custom_handles, fontsize=8, bbox_to_anchor=(1.02, 1), loc='upper left')
                
                plt.grid(True, linestyle=":", alpha=0.7)
                plt.ylim(ymin=0)
                plt.savefig(os.path.join(OUTDIR, f"{arch}_power_model_piecewise.png"), dpi=300, bbox_inches="tight")
                plt.close()
        else:
            print(" -> Skipping: No frequency data found for modeling.")

    print("🎉 Done!")
