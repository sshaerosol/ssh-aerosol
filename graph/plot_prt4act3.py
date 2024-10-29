import os
import numpy as np
import matplotlib.pyplot as plt

# folder_path = 'results/zz/'
folder_path = "../results/prt4act3/"
# Define line styles and colors
line_styles = ["-", "--", "-.", ":"]
colors = ["b", "g", "r", "c", "m", "y", "k"]
# Labels
act3 = ["Ideal", "Unifac", "Unifac+Aiomfac"]
prt4 = ["default", "HPHI", "HPHO", "BOTH"]

filenames = os.listdir(folder_path)
print("\nPlotting for test case prt4act3 ...")

for iact, act in enumerate(act3):
    data = []
    print(f"Reading data and plotting for {act} ...")
    for prt in prt4:
        # Get folder name
        iname = f"{prt}{iact + 1}"
        if iname not in filenames or not os.path.isdir(os.path.join(folder_path, iname)):
            raise ValueError(f"Missing result folder {iname} in {folder_path}")
        # Get data
        ifile = os.path.join(folder_path, iname, "aero", "Organic.txt")
        if not os.path.isfile(ifile):
            raise ValueError(f"Missing data file {ifile}")
        with open(ifile, "r") as f:
            data.append(np.loadtxt(f))

    # Plot
    nt = len(data[0])
    time_steps = np.arange(nt)
    fig, ax = plt.subplots()
    for i, label in enumerate(prt4):
        ax.plot(
            time_steps,
            data[i],
            label=label,
            linestyle=line_styles[i % len(line_styles)],
            color=colors[i % len(colors)],
        )

    # Set limits
    ax.set_xlim(0, nt)
    ax.set_ylim(0, 15)

    ax.set_xlabel("Time Step (#)", fontsize=18)
    ax.set_ylabel(r"Concentration ($\mu$g/m$^3$)", fontsize=18)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="best", framealpha=0.5, fontsize=16)


    ax.set_title(act, fontsize=20)
    plt.savefig(f"{act}.png", dpi=300)
    plt.close()
    print(f"{act}.png created")

# Use montage to combine images
print("Combining images...")
os.system("montage Ideal.png Unifac.png Unifac+Aiomfac.png -geometry +1+1 -tile 3x1 prt4act3.png")
os.system("rm Ideal.png Unifac.png Unifac+Aiomfac.png")
print("prt4act3.png created")
