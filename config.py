"""
Configuration file for spacecraft attitude control simulation.
Contains all physical parameters, limits, and simulation settings.
"""
import numpy as np

# =========================================================
# Spacecraft Physical Parameters
# =========================================================

# Mass parameters
M_TOTAL = 2.6           # total CubeSat mass [kg]
M_RW = 0.13             # reaction wheel mass [kg]

# Geometric parameters
SIDE = 0.1              # CubeSat side length [m]
RW_OFFSET = 0.04        # distance of wheel mass from COM [m]

# Wheel geometric parameters
WHEEL_RADIUS = 0.01     # reaction wheel radius [m]
WHEEL_THICKNESS = 0.005 # reaction wheel thickness [m]

# =========================================================
# Inertia Calculations
# =========================================================

# CubeSat bus inertia (uniform cube)
I_BUS = (1/12) * (M_TOTAL - 4*M_RW) * SIDE**2 * np.eye(3)

# Rotor spin inertia (resists wheel spin)
J_RW = 0.5 * M_RW * (0.02**2)  # rotor spin inertia

# Reaction wheel offset contribution to spacecraft inertia
I_RW_OFFSET = M_RW * RW_OFFSET**2 * np.eye(3)  # offset contribution

# Total spacecraft inertia
I = I_BUS + 4 * I_RW_OFFSET
I_INV = np.linalg.inv(I)

# =========================================================
# Reaction Wheel Configuration
# =========================================================

# Pyramid reaction wheel axes (normalized)
RW_AXES = np.array([[1,1,1], [1,-1,-1], [-1,1,-1], [-1,-1,1]], dtype=float)
RW_AXES = np.array([v/np.linalg.norm(v) for v in RW_AXES])

# =========================================================
# Reaction Wheel Limits
# =========================================================

TAU_MAX = 0.007         # maximum torque per wheel [N·m]
RPM_MAX = 8000.0        # maximum wheel speed [RPM]
W_RW_MAX = RPM_MAX * 2 * np.pi / 60  # maximum wheel speed [rad/s]

# =========================================================
# Controller Parameters
# =========================================================

T_SETTLE = 6.0          # desired settling time [seconds]
ZETA = 0.707            # critical damping ratio

# =========================================================
# Simulation Parameters
# =========================================================

DT = 0.1                # simulation time step [seconds]
T_TOTAL = 10            # total simulation time [seconds]
T_COMMAND = 1.0         # time when attitude command is issued [seconds]
T_FAILURE = 2.0         # time when RW failure occurs [seconds]
FAILURE_DURATION = 4.0  # duration of RW failure [seconds]

# =========================================================
# Desired Attitude (Euler Angles)
# =========================================================

ROLL_DESIRED = np.deg2rad(20)   # desired roll angle [radians]
PITCH_DESIRED = np.deg2rad(90)  # desired pitch angle [radians]
YAW_DESIRED = np.deg2rad(60)    # desired yaw angle [radians]

# =========================================================
# Visualization Parameters
# =========================================================

FIGURE_SIZE = (100, 100)        # figure size for plots
ANIMATION_INTERVAL = 1          # animation interval [ms]
AXIS_LIMIT = 0.15               # 3D axis limits [m]
AXIS_LENGTH = 0.12              # reference frame axis length [m]
