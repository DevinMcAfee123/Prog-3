from metpy.plots import SkewT
from metpy.units import units
import matplotlib.pyplot as plt
import numpy as np
import utils
import metpy.calc as mpcalc

fig_idx =1
for day in [1,4]:
    
    procYear, procMonth, procDay, procHour, procSite = 2011, 7, day , 0, 72357
    soundingInfo, p, h, T, Td, RH, w, wdir, sknt, thta, thte, thtv =         utils.getWyoming(procYear, procMonth, procDay, procHour, procSite)
    p = p[1:]
    h = h[1:]
    T = T[1:]
    Td = Td[1:]
    sknt = sknt[1:]
    wdir = wdir[1:]
    #using poisson's equation to plot parcel temperature at ten evenly spaced pressure levels
    parcel_pressures = np.linspace(p[0], 700, 10)
    parcel_temps = [((T[0]+273.15)*(pressure/p[0])**(2/7))-273.15 for pressure in parcel_pressures]
    parcel_prof = mpcalc.parcel_profile(p*units.hPa, T[0]*units.degC, Td[0]*units.degC).to('degC')
    fig = plt.figure(figsize=(9, 9))
    skew = SkewT(fig, rotation=30)
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    ww = np.array([0.0002, 0.0004, 0.001, 0.002, 0.004, 0.007, 0.01, 0.016, 0.024, 0.032, 0.045, 0.060]).reshape(-1, 1)
    pp = np.linspace(100, 1000) * units.mbar
    skew.plot_mixing_lines(ww, pp)
    skew.plot(p*units.hPa, T*units.degC, 'r', linewidth=1.5)
    skew.plot(p*units.hPa, Td*units.degC, 'g', linewidth=1.5)
    𝚜𝚔𝚎𝚠.𝚙𝚕𝚘𝚝(parcel_pressures*𝚞𝚗𝚒𝚝𝚜.𝚑𝙿𝚊, parcel_temps*𝚞𝚗𝚒𝚝𝚜.𝚍𝚎𝚐𝙲, 'ko', markerfacecolor='black') #plotting discrete parcel temperature values
    skew.plot(p*units.hPa, parcel_prof, 'k', linewidth=2) #plotting parcel temperature profile
    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0]*units.hPa, T[0]*units.degC, Td[0]*units.degC)
    skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')
    plt.xlabel(f'degrees_Celsius\n\nFig {fig_idx} Environmental and Parcel Temperature Profiles')
    plt.title(soundingInfo)
    plt.tight_layout()
    plt.savefig(f'FIG{fig_idx}')

    plt.show()
    fig_idx+=1
    
    p_irr = p[np.where(p>=700)]
    T_irr = [((T[0]+273.15)*(pressure/p[0])**(2/7))-273.15 for pressure in p_irr]
    parcel_pressures = np.linspace(p[0], 700, 10)
    with open(f'FIG{fig_idx}.txt', 'w') as wf:
        for i in range(len(p_irr)):
            wf.write(f'\t\tAt {p_irr[i]:06.2f} hPa, Delta T: {T_irr[i]-T[i]:05.2f} C\n')
        wf.write(f'\n\tFig {fig_idx} Parcel Temperature - Environmental Temperature\n')
    fig_idx+=1
    
    thta = thta[1:]
    fig = plt.figure(figsize=(9, 9))
    skew = SkewT(fig, rotation=30)
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    ww = np.array([0.0002, 0.0004, 0.001, 0.002, 0.004, 0.007, 0.01, 0.016, 0.024, 0.032, 0.045, 0.060]).reshape(-1, 1)
    pp = np.linspace(100, 1000) * units.mbar
    skew.plot_mixing_lines(ww, pp)
    skew.plot(p_irr, thta[:len(p_irr)]-273.15, 'k', linewidth=1.5)
    plt.xlabel(f'degrees_Celsius\n\nFig {fig_idx} Potential Temperature from Surface to 700 hPa')
    plt.tight_layout()
    plt.title(soundingInfo)
    plt.savefig(f'FIG{fig_idx}')
    plt.show()
    fig_idx+=1

