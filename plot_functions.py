"""
Plotting and animation functions for spacecraft visualization.
"""
import numpy as np
import matplotlib.pyplot as plt
from math_functions import quat_to_dcm
from config import FIGURE_SIZE, AXIS_LIMIT, AXIS_LENGTH


def plot_rw(ax, pos, axis, radius=0.01, thickness=0.005, color='purple'):
    """
    Plot a reaction wheel as a cylinder along `axis` at position `pos`.
    
    Parameters:
    -----------
    ax : matplotlib 3D axis
        The axis to plot on
    pos : array
        Position of the wheel center
    axis : array
        Axis direction of the wheel
    radius : float
        Wheel radius
    thickness : float
        Wheel thickness
    color : str
        Color of the wheel
    """
    # Number of points for circle
    n = 20
    theta = np.linspace(0, 2*np.pi, n)
    circle = np.array([radius*np.cos(theta), radius*np.sin(theta), np.zeros(n)])  # z=0
    
    # Build top and bottom circles
    top = circle + np.array([[0],[0],[thickness/2]])
    bottom = circle - np.array([[0],[0],[thickness/2]])
    
    # Align cylinder with axis
    z_axis = np.array([0,0,1])
    axis = axis / np.linalg.norm(axis)
    v = np.cross(z_axis, axis)
    s = np.linalg.norm(v)
    c = np.dot(z_axis, axis)
    if s < 1e-6:
        rot = np.eye(3)
    else:
        vx = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
        rot = np.eye(3) + vx + vx@vx*(1/(1+c))
    
    top_rot = (rot @ top).T + pos
    bottom_rot = (rot @ bottom).T + pos
    
    # Connect top and bottom
    for i in range(n):
        xs = [top_rot[i,0], bottom_rot[i,0]]
        ys = [top_rot[i,1], bottom_rot[i,1]]
        zs = [top_rot[i,2], bottom_rot[i,2]]
        ax.plot(xs, ys, zs, color=color)


def plot_rw_pyramid(ax, w_rw, rw_axes, scale=0.005):
    """
    Plot reaction wheels as vectors along pyramid axes, length scaled by speed.
    
    Parameters:
    -----------
    ax : matplotlib 3D axis
        The axis to plot on
    w_rw : array
        Wheel speeds
    rw_axes : array
        Reaction wheel axes (4x3)
    scale : float
        Scaling factor for visualization
    """
    for i, axis in enumerate(rw_axes):
        vec = axis / np.linalg.norm(axis) * scale * w_rw[i]
        color = 'red' if w_rw[i] >= 0 else 'blue'
        ax.plot([0, vec[0]], [0, vec[1]], [0, vec[2]], color=color, lw=3)


def create_animation_update_function(
    t, state_hist, qd_hist, failure_hist, w_rw_hist, tau_rw_hist, trace_err_hist,
    wd_hist, side, rw_axes, rw_offset, w_rw_max, tau_max
):
    """
    Create the update function for the animation.
    
    Parameters:
    -----------
    t : array
        Time array
    state_hist : array
        State history
    qd_hist : array
        Desired quaternion history
    failure_hist : array
        Failure history
    w_rw_hist : array
        Wheel speed history
    tau_rw_hist : array
        Wheel torque history
    trace_err_hist : array
        DCM trace error history
    wd_hist : array
        Desired angular rate history
    side : float
        CubeSat side length
    rw_axes : array
        Reaction wheel axes
    rw_offset : float
        Reaction wheel offset from center
    w_rw_max : float
        Maximum wheel speed
    tau_max : float
        Maximum torque
    
    Returns:
    --------
    function : Update function for FuncAnimation
    """
    # Cube geometry
    cube_pts = np.array([
        [-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
        [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]
    ]) * side / 2
    cube_edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
    
    # Pyramid edges
    edges = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    
    def update(k, ax_anim, ax_trace, ax_w, ax_rw, ax_tau):
        """Animation update function."""
        ax_anim.cla()
        ax_anim.set_xlim([-AXIS_LIMIT, AXIS_LIMIT])
        ax_anim.set_ylim([-AXIS_LIMIT, AXIS_LIMIT])
        ax_anim.set_zlim([-AXIS_LIMIT, AXIS_LIMIT])
        ax_anim.set_box_aspect([1,1,1])
        ax_anim.set_title(f"CubeSat Attitude   t={t[k]:.1f}s    RW1 Failed: {failure_hist[k]}")
        
        C = quat_to_dcm(state_hist[k,:4])
        Cd = quat_to_dcm(qd_hist[k])
        
        # Cube edges
        pts = (C @ cube_pts.T).T
        for e in cube_edges:
            ax_anim.plot(*zip(pts[e[0]], pts[e[1]]), color='k')
        
        # Reaction wheel positions
        rw_positions = rw_axes * rw_offset
        # Rotate wheel positions with current CubeSat orientation
        rw_rot = (C @ rw_positions.T).T  # shape (4,3)

        # Plot pyramid satellite as simple lines
        for i in range(4):
            p = C @ rw_axes[i] * rw_offset
            ax_anim.scatter(p[0], p[1], p[2], color='purple', s=20)
        
        # Connect all edges
        for e in edges:
            ax_anim.plot([rw_rot[e[0],0], rw_rot[e[1],0]],
                         [rw_rot[e[0],1], rw_rot[e[1],1]],
                         [rw_rot[e[0],2], rw_rot[e[1],2]], color='purple')
        
        # Reference and current frames
        for i, c in enumerate(['r','g','b']):
            ax_anim.plot([0, Cd[0,i]*AXIS_LENGTH], [0, Cd[1,i]*AXIS_LENGTH], [0, Cd[2,i]*AXIS_LENGTH],
                         color=c, alpha=0.3, lw=2)
            ax_anim.plot([0, C[0,i]*AXIS_LENGTH], [0, C[1,i]*AXIS_LENGTH], [0, C[2,i]*AXIS_LENGTH],
                         color=c, lw=3)
        
        # DCM trace error
        ax_trace.cla()
        ax_trace.plot(t[:k], trace_err_hist[:k], label='trace(C_err)')
        ax_trace.set_title("DCM Error Trace")
        ax_trace.legend()
        
        # Angular rates
        ax_w.cla()
        ax_w.plot(t[:k], state_hist[:k,4:], label=['wx','wy','wz'])
        ax_w.plot(t[:k], wd_hist[:k], '--', label='desired')
        ax_w.set_title("Angular Rates")
        ax_w.legend()
        ax_w.set_ylabel("rad/s")
        
        # Reaction wheel speeds
        ax_rw.cla()
        ax_rw.plot(t[:k], w_rw_hist[:k], label=[f'RW{i+1}' for i in range(4)])
        ax_rw.axhline(w_rw_max, color='k', ls='--')
        ax_rw.axhline(-w_rw_max, color='k', ls='--')
        ax_rw.set_title("Reaction Wheel Speeds")
        ax_rw.set_ylim([-80,80])
        ax_rw.legend()
        ax_rw.set_ylabel("rad/s")
        
        # Reaction wheel torques
        ax_tau.cla()
        ax_tau.plot(t[:k], tau_rw_hist[:k], label=[f'RW{i+1}' for i in range(4)])
        ax_tau.axhline(tau_max, color='k', ls='--')
        ax_tau.axhline(-tau_max, color='k', ls='--')
        ax_tau.set_title("Reaction Wheel Torques")
        ax_tau.legend()
        ax_tau.set_ylabel("N.m")
    
    return update


def setup_figure():
    """
    Set up the figure with subplots for animation.
    
    Returns:
    --------
    tuple : (fig, ax_anim, ax_trace, ax_w, ax_rw, ax_tau)
    """
    fig = plt.figure(figsize=FIGURE_SIZE)
    gs = fig.add_gridspec(2, 3, width_ratios=[2,1,1], height_ratios=[1,1], figure=fig)
    
    ax_anim = fig.add_subplot(gs[:,0], projection='3d')
    ax_trace = fig.add_subplot(gs[0,1])
    ax_w = fig.add_subplot(gs[0,2])
    ax_rw = fig.add_subplot(gs[1,1])
    ax_tau = fig.add_subplot(gs[1,2])
    
    return fig, ax_anim, ax_trace, ax_w, ax_rw, ax_tau
