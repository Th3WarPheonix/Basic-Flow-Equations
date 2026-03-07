import international_atmosphere as intlatm
import numpy as np
import matplotlib.pyplot as plt
import pint

ureg = pint.UnitRegistry()

def get_std_day(day, units):
    LOWER_ALT = 0 
    HIGHER_ALT = 80000
    NUM_POINTS = 1000 # Number of altitudes at which to plot the data

    altitude = np.linspace(LOWER_ALT, HIGHER_ALT, NUM_POINTS) * ureg.meters

    if day == 'standard':
        temperatures = intlatm.StandardDay(altitude.magnitude).temperature * ureg.kelvin
        pressures = intlatm.StandardDay(altitude.magnitude).pressure * ureg.pascal
        densities = intlatm.StandardDay(altitude.magnitude).density * ureg.kilogram / ureg.meter**3
    elif day == 'hot':
        temperatures = intlatm.HotDay(altitude.magnitude).temperature * ureg.kelvin
        pressures = intlatm.HotDay(altitude.magnitude).pressure * ureg.pascal
        densities = intlatm.HotDay(altitude.magnitude).density * ureg.kilogram / ureg.meter**3
    elif day == 'cold':
        temperatures = intlatm.ColdDay(altitude.magnitude).temperature * ureg.kelvin
        pressures = intlatm.ColdDay(altitude.magnitude).pressure * ureg.pascal
        densities = intlatm.ColdDay(altitude.magnitude).density * ureg.kilogram / ureg.meter**3
    elif day == 'tropical':
        temperatures = intlatm.TropicalDay(altitude.magnitude).temperature * ureg.kelvin
        pressures = intlatm.TropicalDay(altitude.magnitude).pressure * ureg.pascal
        densities = intlatm.TropicalDay(altitude.magnitude).density * ureg.kilogram / ureg.meter**3

    if units == 'metric':
        return altitude, temperatures, pressures, densities
    else:
        altitude.ito(ureg.feet)
        temperatures.ito(ureg.fahrenheit)
        pressures.ito(ureg.psi)
        densities.ito(ureg.pound/ureg.ft**3)
        return altitude, temperatures, pressures, densities

def plot_properties(altitude, temperature, pressure, density, units):

    if units == 'metric':
        Y_LABEL = 'Altitude (m)'
        TEMP_LABEL = 'Temperature ($K$)'
        PRESSURE_LABEL = 'Pressure ($N/m^2$)'
        DENSITY_LABEL = 'Density ($kg/m^3$)'
    else:
        Y_LABEL = 'Altitude (ft)'
        TEMP_LABEL = 'Temperature ($F$)'
        PRESSURE_LABEL = 'Pressure ($psi$)'
        DENSITY_LABEL = 'Density ($lbm/ft^3$)'

    if day == 'standard':
        FIG_FILENAME = 'standard_day_atmosphere.png'
        FIG_TITLE = 'Standard Day'
    elif day == 'hot':
        FIG_FILENAME = 'hot_day_atmosphere.png'
        FIG_TITLE = 'Hot Day'
    elif day == 'cold':
        FIG_FILENAME = 'cold_day_atmosphere.png'
        FIG_TITLE = 'Cold Day'
    elif day == 'tropical':
        FIG_FILENAME = 'tropical_day_atmosphere.png'
        FIG_TITLE = 'Tropical Day'

    AX1_COLOR = 'red'
    AX2_COLOR = 'blue'
    AX3_COLOR = 'darkgreen'
    PLOT_PADDING_BOTTOM = .15
    PLOT_PADDING_TOP = .95
    SPACING = -PLOT_PADDING_BOTTOM/3-.02
    X_MARGINS = .01 # lines up start and end of x ticks

    fig, ax1 = plt.subplots(figsize=(9, 7))
    plt.subplots_adjust(top=PLOT_PADDING_TOP, bottom=PLOT_PADDING_BOTTOM)

    ax1.plot(temperature, altitude, label=TEMP_LABEL, color=AX1_COLOR, zorder=15)
    ax1.set_ylabel(Y_LABEL)
    temp_xticks = np.linspace(min(temperature.magnitude), max(temperature.magnitude), 10)
    ax1.set_xticks(temp_xticks)
    ax1.set_xticklabels(np.around(temp_xticks, 0))
    ax1.tick_params(axis='x', colors=AX1_COLOR, zorder=0)
    ax1.margins(x=X_MARGINS)

    ax2 = ax1.twiny()
    ax2.plot(pressure, altitude, label=PRESSURE_LABEL, color=AX2_COLOR, zorder=15)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.spines['bottom'].set_position(("axes", SPACING))
    ax2.set_xticks(np.around(np.linspace(0, max(pressure.magnitude), 10), 0))
    ax2.tick_params(axis='x', colors=AX2_COLOR, zorder=0)
    ax2.margins(x=X_MARGINS)

    ax3 = ax1.twiny()
    ax3.plot(density, altitude, label=DENSITY_LABEL, color=AX3_COLOR, zorder=15)
    ax3.xaxis.set_ticks_position('bottom')
    ax3.spines['bottom'].set_position(("axes", 2*SPACING))
    ax3.set_xticks(np.around(np.linspace(0, max(density.magnitude), 10), 3))
    ax3.tick_params(axis='x', colors=AX3_COLOR, zorder=0)
    ax3.margins(x=X_MARGINS)

    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    handles3, labels3 = ax3.get_legend_handles_labels()

    fig.legend(handles1+handles2+handles3, labels1+labels2+labels3, framealpha=1)

    ax1.grid(axis='both')
    ax1.minorticks_on()
    ax2.minorticks_on()
    ax3.minorticks_on()
    ax1.grid(which='minor', linestyle=':')
    ax1.set_ylim(0)
    fig.suptitle(FIG_TITLE)
    plt.savefig(r'./reference_figures/' + FIG_FILENAME)

if __name__ == '__main__':
    day = 'cold'
    units = 'english'
    atpd = get_std_day(day, units)
    plot_properties(atpd[0], atpd[1], atpd[2], atpd[3], units)
