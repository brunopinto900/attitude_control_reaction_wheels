"""
Refactored spacecraft attitude control simulation.
Main script that orchestrates the simulation using modular components.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib
matplotlib.use('TkAgg')

# Import configuration
from config import*

# Import mathematical functions
from math_functions import (
    euler_to_quat, quat_mult, quat_conj, quat_norm, quat_to_dcm
)

# Import controller
from controller import cascaded_attitude_controller

# Import dynamics functions
from dynamics import rk4

# Import reaction wheel functions
from reaction_wheels import addAlignmentErrors, update_reaction_wheels

# Import plot functions
from plot_functions import setup_figure, create_animation_update_function


# =========================================================
# Derived Configuration
# =========================================================

# Reaction wheel axis matrices
RW_axis = RW_AXES.T
RW_axis_pinv = np.linalg.pinv(RW_axis)
RW_axis_with_err = addAlignmentErrors(RW_AXES).T

# Compute wheel center positions
rw_pos = RW_AXES * RW_OFFSET


# =========================================================
# Main Simulation Function
# =========================================================

def main():
    """Main simulation function with refactored code."""
    # Simulation parameters
    N = int(T_TOTAL / DT)
    t = np.arange(N) * DT

    # Initialize state
    state = np.zeros(7)
    state[0] = 1.0  # initial quaternion
    w_rw = np.zeros(4)

    # Desired attitude quaternion
    qd_cmd = euler_to_quat(ROLL_DESIRED, PITCH_DESIRED, YAW_DESIRED)
    qd_cmd = qd_cmd / np.linalg.norm(qd_cmd)

    # History arrays
    failure_hist = []
    state_hist = []
    qd_hist = []
    wd_hist = []
    w_rw_hist = []
    tau_rw_hist = []
    trace_err_hist = []

    # Simulation loop
    for k in range(N):
        t_now = t[k]

        # Set desired state
        if t[k] < T_COMMAND:
            qd, wd = np.array([1,0,0,0]), np.zeros(3)
        else:
            qd, wd = qd_cmd, np.zeros(3)
        
        # Extract current state
        q, w = state[:4], state[4:]
        q_err = quat_mult(quat_conj(qd), q)
        if q_err[0] < 0:
            q_err = -q_err
        C_err = quat_to_dcm(q_err)
        trace_err_hist.append(np.trace(C_err))
        
        # Compute control torque using cascaded controller
        torque_cmd, wd = cascaded_attitude_controller(
            q_err[1:], w, I,
            Ts=T_SETTLE,
            zeta=ZETA,
            omega_max=2.0,
            torque_max=None
        )
        
        # Map to wheel torques
        tau_rw_cmd = np.clip(RW_axis_pinv @ torque_cmd, -TAU_MAX, TAU_MAX)
        
        # Check for failure
        failure = (t_now >= T_FAILURE) and (t_now <= (T_FAILURE + FAILURE_DURATION))

        # Update reaction wheels
        w_rw, tau_rw_actual = update_reaction_wheels(
            w_rw, tau_rw_cmd, J_RW, W_RW_MAX, DT, failure
        )
        
        # Integrate spacecraft dynamics
        state = rk4(state, RW_axis_with_err @ tau_rw_actual, DT, I_INV, I)
        state[:4] = quat_norm(state[:4])
        
        # Store history
        failure_hist.append(failure)
        state_hist.append(state.copy())
        qd_hist.append(qd.copy())
        wd_hist.append(wd.copy())
        w_rw_hist.append(w_rw.copy())
        tau_rw_hist.append(tau_rw_actual.copy())

    # Convert to arrays
    failure_hist = np.array(failure_hist)
    state_hist = np.array(state_hist)
    qd_hist = np.array(qd_hist)
    wd_hist = np.array(wd_hist)
    w_rw_hist = np.array(w_rw_hist)
    tau_rw_hist = np.array(tau_rw_hist)
    trace_err_hist = np.array(trace_err_hist)

    # =========================================================
    # Visualization
    # =========================================================

    fig, ax_anim, ax_trace, ax_w, ax_rw, ax_tau = setup_figure()

    # Create update function for animation
    update_func = create_animation_update_function(
        t, state_hist, qd_hist, failure_hist, w_rw_hist, tau_rw_hist,
        trace_err_hist, wd_hist, SIDE, RW_AXES, RW_OFFSET, W_RW_MAX, TAU_MAX
    )

    # Animation wrapper
    def animate(k):
        update_func(k, ax_anim, ax_trace, ax_w, ax_rw, ax_tau)

    ani = FuncAnimation(fig, animate, frames=N, interval=ANIMATION_INTERVAL)

    print("Simulation complete. Displaying animation...")
    plt.show()


if __name__ == "__main__":
    main()

