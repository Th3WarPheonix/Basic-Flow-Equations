from ambiance import Atmosphere
import numpy as np
import matplotlib.pyplot as plt


def plotProps():
    plt.style.use('bmh')

    LOWER_ALT = 0 # Lowest altitude to graph
    HIGHER_ALT = 80000 # Highest altitude to graph
    NUM_POINTS = 1000 # Number of altitudes at which to plot the data

    heights = np.linspace(LOWER_ALT, HIGHER_ALT, NUM_POINTS)

    temperatures = Atmosphere(heights).temperature
    pressures = Atmosphere(heights).pressure
    densities = Atmosphere(heights).density

    Y_LABEL = 'Altitude (m)'
    AX1_COLOR = 'red'
    AX2_COLOR = 'blue'
    AX3_COLOR = 'darkgreen'
    PLOT_PADDING_BOTTOM = .15
    PLOT_PADDING_TOP = .95
    SPACING = -PLOT_PADDING_BOTTOM/3-.02
    X_MARGINS = .01

    fig=plt.figure(figsize=(9, 7))
    plt.subplots_adjust(top=PLOT_PADDING_TOP, bottom=PLOT_PADDING_BOTTOM)
    ax1=fig.add_subplot(111, label="1")
    ax2=fig.add_subplot(111, label="2", frame_on=False)
    ax3=fig.add_subplot(111, label="3", frame_on=False)

    ax1.plot(temperatures, heights, label='Temperature ($K$)', color=AX1_COLOR, zorder=15)
    ax1.set_ylabel(Y_LABEL)
    #ax1.set_xlabel("Temperatures ($K$)", color=AX1_COLOR, loc='left')
    ax1.set_xticks(np.linspace(min(temperatures), max(temperatures), 10))
    # Allows user to arbitrarily set axis tick values, therefore they show up rounded when in reality they are not (comment out ths line to see what happens)
    ax1.set_xticklabels(np.around(np.linspace(min(temperatures), max(temperatures), 10), 0))
    ax1.tick_params(axis='x', colors=AX1_COLOR, zorder=0)
    ax1.margins(x=X_MARGINS)

    ax2.plot(pressures, heights, label='Pressure ($N/m^2$)', color=AX2_COLOR, zorder=15)
    ax2.spines.bottom.set_position(("axes", SPACING))
    #ax2.set_xlabel('Pressure ($N/m^2$)', color=AX2_COLOR, loc='left')
    ax2.set_xticks(np.around(np.linspace(0, max(pressures), 10), 0))
    ax2.tick_params(axis='x', colors=AX2_COLOR, zorder=0)
    ax2.margins(x=X_MARGINS)
    ax2.set_yticks([])

    ax3.plot(densities, heights, label='Density ($kg/m^3$)', color=AX3_COLOR, zorder=15)
    ax3.spines.bottom.set_position(("axes", 2*SPACING))
    #ax3.set_xlabel('Density ($kg/m^3$)', color=AX3_COLOR, loc='left')
    ax3.set_xticks(np.around(np.linspace(0, max(densities), 10), 3))
    ax3.tick_params(axis='x', colors=AX3_COLOR, zorder=0)
    ax3.margins(x=X_MARGINS)
    ax3.set_yticks([])

    plot_handles = []
    plot_labels = []
    handles, labels = ax1.get_legend_handles_labels()
    plot_handles.append(handles[0])
    plot_labels.append(labels[0])
    handles, labels = ax2.get_legend_handles_labels()
    plot_handles.append(handles[0])
    plot_labels.append(labels[0])
    handles, labels = ax3.get_legend_handles_labels()
    plot_handles.append(handles[0])
    plot_labels.append(labels[0])
    fig.legend(plot_handles, plot_labels, loc='right', framealpha=1)

    #ax1.grid(alpha=.5, zorder=0)

    plt.title("Properties Through Earth's Atmoshpere")
    #plt.savefig('Atmospheric Properties')
    plt.show()

if __name__ == '__main__':
    plotProps()
    """
    https://pypi.org/project/ambiance/
    List of properties
    Collision frequency (collision_frequency)
    Density (density)
    Dynamic viscosity (dynamic_viscosity)
    Geometric height above MSL (h)
    Geopotential height (H)
    Gravitational acceleration (grav_accel)
    Kinematic viscosity (kinematic_viscosity)
    Layer names (layer_name) [string array]
    Mean free path (mean_free_path)
    Mean particle speed (mean_particle_speed)
    Number density (number_density)
    Pressure (pressure)
    Pressure scale height (pressure_scale_height)
    Specific weight (specific_weight)
    Speed of sound (speed_of_sound)
    Temperature (temperature, temperature_in_celsius)
    Thermal conductivity (thermal_conductivity)
    """