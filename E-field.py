import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation







def elliptical_polarization(t, A, theta, delta_phi):
    omega = 2 * np.pi  # Angular frequency
    Ax = A * np.cos(theta)
    Ay = A * np.sin(theta)
    Ex = Ax * np.cos(omega * t)
    Ey = Ay * np.cos(omega * t + delta_phi)
    return Ex, Ey



def calculate_psi(A, theta, phi):
    numerator = 2 * A**2 * np.cos(theta) * np.sin(theta) * np.cos(phi)
    denominator = A**2 * np.cos(theta)**2 - A**2 * np.sin(theta)**2
    psi = 0.5 * np.arctan(numerator / denominator)
    psi = np.rad2deg(psi)
    return psi





A = 1.0  # Amplitude
theta = np.pi / 4
phi = np.pi / 4  # +-, Phase difference between x and y component. pi/2 for circle


psi = round(calculate_psi(A, theta, phi),1)






def update(frame):
    t = frame / 30.0  # Adjust for smoother animation
    Ex, Ey = elliptical_polarization(t, A, theta, phi)
    trajectory_x.append(Ex)
    trajectory_y.append(Ey)

    plt.clf()
    plt.quiver(0, 0, Ex, Ey, angles='xy', scale_units='xy', scale=1, color='b', width=0.01)
    plt.plot(trajectory_x, trajectory_y, color='r', linestyle='dashed', linewidth=2)

    plt.xlim(-A, A)
    plt.ylim(-A, A)
    plt.xlabel('$E_x$')
    plt.ylabel('$E_y$')
    
    plt.title(f'Time: {t:.1f} seconds. Angle ψ: {psi} deg')


trajectory_x = []
trajectory_y = []

fig = plt.figure()
ani = FuncAnimation(fig, update, frames=600, interval=20, blit=False) 
plt.show()
