# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 02:54:35 2025

@author: ASUS
"""

import numpy as np
import matplotlib.pyplot as plt

# Physical / numerical parameters
dt = 1e-18
a = 2.4*1e-10            # lattice constant
dx = dy = dz = a/20
c = 3e8
gamma = 0.3 * ((2 * np.pi * c)/a)
sigma = 10 * ((c**2 * 4 * np.pi**2) / (a**2))
N = 1000               # fewer steps for demonstration
C0 = 1
deltaV = dx * dy * dz
omega_0 = 0
mu = 1.25663706e-6
Ks = np.sqrt(((48*np.pi**2*C0**2)/(N*deltaV)) * sigma * gamma)
epsilon_inf = 1

# Grid definition
nx, ny, nz = 40, 40, 20
Lx = nx * dx
Ly = ny * dy
Lz = nz * dz
x = np.linspace(0, Lx, nx, endpoint=False)
y = np.linspace(0, Ly, ny, endpoint=False)
z = np.linspace(0, Lz, nz, endpoint=False)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Curl on grid
def curl(F, dx=dx, dy=dy, dz=dz):
    Fx, Fy, Fz = F[...,0], F[...,1], F[...,2]
    dFz_dy = np.gradient(Fz, dy, axis=1)
    dFy_dz = np.gradient(Fy, dz, axis=2)
    curl_x = dFz_dy - dFy_dz

    dFx_dz = np.gradient(Fx, dz, axis=2)
    dFz_dx = np.gradient(Fz, dx, axis=0)
    curl_y = dFx_dz - dFz_dx

    dFy_dx = np.gradient(Fy, dx, axis=0)
    dFx_dy = np.gradient(Fx, dy, axis=1)
    curl_z = dFy_dx - dFx_dy

    return np.stack((curl_x, curl_y, curl_z), axis=-1)

# Simulation on grid
def sim_grid(P0, E0, H0):
    sol = []

    # Initialize
    P_nm1 = P0.copy()
    P_n   = P0.copy()
    E_n   = E0.copy()
    H_nm1 = H0.copy()
    H_n = H0.copy()

    denom = 2 + gamma * dt

    for i in range(N):
        # random K field
        K_t = np.random.uniform(-Ks/2,Ks/2, size=P0.shape)

        # P update
        P_np1 = (2 * dt**2 / denom) * (
            sigma * E_n + K_t
            + (2/dt**2 - omega_0**2) * P_n
            + (gamma/(2*dt) - 1/dt**2) * P_nm1
        )

        # E update
        E_np1 = E_n + (4 * np.pi / epsilon_inf) * (P_n - P_np1)

        # H update
        curl_E = curl(E_np1)
        H_np1 = H_nm1 - (2 * dt / mu) * curl_E

        # Poynting vector
        S = np.cross(E_np1, H_np1)

        sol.append(S)

        # Shift for next step
        P_nm1, P_n = P_n, P_np1
        E_n         = E_np1
        H_nm1,H_n       = H_n, H_np1

    return np.array(sol)

# Initial conditions: P and H zero, E with spatial variation
P0 = np.zeros((nx, ny, nz, 3))
H0 = np.zeros((nx, ny, nz, 3))
E0 = np.zeros((nx, ny, nz, 3))
# Run simulation
result = sim_grid(P0, E0, H0)
result = result[1000::,:,:,:,:]

# Select a mid‐plane slice for visualization
k = 0  #z-slice
S_slice = result[97,:, :, k, :]  # shape (nx, ny, 3)
X_slice = X[:, :, k]
Y_slice = Y[:, :, k]

Sx = S_slice[:, :, 0]
Sy = S_slice[:, :, 1]
Sz = S_slice[:,:,2]

# magnitude
S_mag = np.sqrt(Sx**2 + Sy**2 + Sz**2)

# contour plot
plt.figure(figsize=(6,5))
cs = plt.contourf(X_slice, Y_slice, S_mag)  # default colormap
plt.colorbar(cs, label='|S| (units)')
plt.title(f'Contour of Poynting Vector Magnitude |S| (z={k*dz})')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.tight_layout()
plt.show()

# === Reconstruct your grid from dx,dy,dz and the shape of result ===
# result.shape == (time_steps, nx, ny, nz, 3)
time_steps, nx, ny, nz, _ = result.shape
x = np.arange(nx) * dx
y = np.arange(ny) * dy
z = np.arange(nz) * dz
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# === Parameters ===
t_index = 99    # which time slice to visualize
step    = 4    # subsample every 4th point for clarity

# === Subsample the grid and Poynting vectors ===
X3 = X[::step, ::step, ::step]
Y3 = Y[::step, ::step, ::step]
Z3 = Z[::step, ::step, ::step]
S3 = result[t_index, ::step, ::step, ::step, :]  # shape (nx', ny', nz', 3)
U, V, W = S3[...,0], S3[...,1], S3[...,2]

# === 3D quiver plot ===
fig = plt.figure(figsize=(8, 6))
ax  = fig.add_subplot(111, projection='3d')
ax.quiver(
    X3, Y3, Z3,          # vector origins
    U, V, W,             # vector components
    length=dx*2,         # scale arrow length
    normalize=True,      # normalize arrows for uniform length
    linewidth=0.5
)

ax.set_title('3D Poynting Vector Field (Subsampled)')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('z (m)')
plt.tight_layout()
plt.show()

sum_over_xyz = result.sum(axis=(1,2,3))
Flux_t = np.einsum('ij,j',sum_over_xyz,dx*dy*np.array([0,0,1]))

fft_vals  = np.fft.fft(Flux_t)
fft_freqs = np.fft.fftfreq(len(Flux_t), dt)
# 3) Plot the magnitude spectrum
plt.figure()
plt.plot(np.abs(fft_freqs/(c/a)), np.abs(fft_vals))
plt.xlabel('Frequency (c/a)')
plt.xlim(0,30)
plt.ylabel('Magnitude |P|')
plt.title('A')
plt.show()