import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


fatras = uproot.open("debug_fatras/hits.root")["hits"].arrays(library="pandas")
geant4 = uproot.open("debug_geant4/hits.root")["hits"].arrays(library="pandas")

fatras.sort_values(["event_id", "index"], inplace=True)
geant4.sort_values(["event_id", "index"], inplace=True)

print(fatras)
print(geant4)

merge = pd.merge(fatras.add_prefix("fatras_"), geant4.add_prefix("geant4_"), how="outer", left_on=["fatras_event_id", "fatras_index"], right_on=["geant4_event_id", "geant4_index"])

print(merge[["fatras_tx", "geant4_tx"]])
print(merge[["fatras_volume_id", "geant4_volume_id"]])
print(merge[["fatras_layer_id", "geant4_layer_id"]])
print(merge[["fatras_sensitive_id", "geant4_sensitive_id"]])

fig = plt.figure("position")
axes = fig.subplots(2,2)

axes[0,0].plot(merge["fatras_tz"], merge["fatras_ty"], "o", linestyle="", label="fatras")
axes[0,0].plot(merge["geant4_tz"], merge["geant4_ty"], "o", linestyle="", label="geant4")
axes[0,0].set_xlabel("z [mm]")
axes[0,0].set_ylabel("y [mm]")
axes[0,0].legend()

axes[0,1].plot(merge["fatras_tx"], merge["fatras_ty"], "o", linestyle="", label="fatras")
axes[0,1].plot(merge["geant4_tx"], merge["geant4_ty"], "o", linestyle="", label="geant4")
axes[0,1].set_xlabel("x [mm]")
axes[0,1].set_ylabel("y [mm]")
axes[0,1].legend()

axes[1,0].plot(merge["fatras_tz"], merge["fatras_tx"], "o", linestyle="", label="fatras")
axes[1,0].plot(merge["geant4_tz"], merge["geant4_tx"], "o", linestyle="", label="geant4")
axes[1,0].set_xlabel("z [mm]")
axes[1,0].set_ylabel("x [mm]")
axes[1,0].legend()

"""
axes[1,1].plot(merge["fatras_index"], merge["fatras_tx"]-merge["geant4_tx"], "o", linestyle="", label="x")
axes[1,1].plot(merge["fatras_index"], merge["fatras_ty"]-merge["geant4_ty"], "o", linestyle="", label="y")
axes[1,1].plot(merge["fatras_index"], merge["fatras_tz"]-merge["geant4_tz"], "o", linestyle="", label="z")
axes[1,1].set_xlabel("fatras - geant4 distance [mm]")
axes[1,1].set_ylabel("hit index")
axes[1,1].legend()
"""

fig = plt.figure("momentum")
axes = fig.subplots(1,3)

fatras_r = np.hypot(merge["fatras_tx"], merge["fatras_ty"], merge["fatras_tz"])
geant4_r = np.hypot(merge["geant4_tx"], merge["geant4_ty"], merge["geant4_tz"])

axes[0].plot(fatras_r, merge["fatras_tpx"], "o", linestyle="", label="fatras")
axes[0].plot(geant4_r, merge["geant4_tpx"], "o", linestyle="", label="geant4")
axes[0].legend()

axes[1].plot(fatras_r, merge["fatras_tpy"], "o", linestyle="", label="fatras")
axes[1].plot(geant4_r, merge["geant4_tpy"], "o", linestyle="", label="geant4")
axes[1].legend()

axes[2].plot(fatras_r, merge["fatras_tpz"], "o", linestyle="", label="fatras")
axes[2].plot(geant4_r, merge["geant4_tpz"], "o", linestyle="", label="geant4")
axes[2].legend()

plt.show()
