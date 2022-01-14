# 530 Stellar Atsmopheres HW 1

import sys
import os
import astropy
from astropy.modeling import models
from astropy import units as u
from astropy import constants as const
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statmospheres import Planck_array_at_temp

def main():
    num_points=100
    wave_numbers=np.linspace(0,12,num_points)/u.um
    input=(wave_numbers*const.c).to(u.Hz)
    Ts=[10000,7000,3000]*u.K
    df=pd.DataFrame({'wave_numbers':wave_numbers})
    for plot in range(3):
        for T in Ts:
            my_B_nus = Planck_array_at_temp(T=T,nu_array=input)
            my_B_nus = [B.value for B in my_B_nus]
            df["my_B_nus"]=my_B_nus
            df=df[~np.isnan(df.my_B_nus)].reset_index(drop=True)
            wave_numbers = [w.value for w in wave_numbers]
            plt(df.wave_numbers,df.my_B_nus)
        if plot == 0:
            x = y = 'lin'
        if plot == 1:
            x = 'log'
            plt.xscale('log')
            y = 'lin'
        if plot == 2:
            x = y = 'log'
            plt.xscale('log')
            plt.yscale('log')
        path_png=os.getcwd()+f'/HW_01_plots/{x}_{y}.png'
        plt.savefig(path_png, bbox_inches='tight')
    return None
# run it!
if __name__ == "__main__":
    main()