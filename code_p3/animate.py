import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import cm


def parse_file(filename):
    with open(filename, "r") as file:
        num_particles, num_timesteps = map(int, file.readline().split())
        timesteps = []
        positions = []
        while True:
            timestep_line = file.readline().strip()
            if not timestep_line:
                break
            timesteps.append(float(timestep_line))
            timestep_positions = []
            for _ in range(num_particles):
                x, y, z, vx, vy, vz = map(float, file.readline().split())
                timestep_positions.append([x, y, z])
            positions.append(np.array(timestep_positions))
    return np.array(timesteps), np.array(positions), num_particles


def animate_particles(
    filename, first_timestep, last_timestep, fps=30, frame_skip=1, zoom=500
):
    timesteps, positions, num_particles = parse_file(filename)
    selected_positions = positions[first_timestep:last_timestep]
    selected_timesteps = timesteps[first_timestep:last_timestep]

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    # zooming in on the origin
    ax.set_xlim(-zoom, zoom)
    ax.set_ylim(-zoom, zoom)
    ax.set_zlim(-zoom, zoom)

    colors = cm.rainbow(np.linspace(0, 1, num_particles))  # Colormap (rainbow)

    scat = ax.scatter([], [], [], s=10)

    interval = 1000 / fps  # in milliseconds

    def update(frame):
        # every Nth frame
        actual_frame = frame * frame_skip

        if actual_frame >= len(selected_positions):
            print(actual_frame, "is out of bounds btw")
            return (scat,)

        scat._offsets3d = (
            selected_positions[actual_frame][:, 0],
            selected_positions[actual_frame][:, 1],
            selected_positions[actual_frame][:, 2],
        )
        scat.set_color(colors)

        timestep_value = selected_timesteps[actual_frame]
        if timestep_value.is_integer():
            timestep_value = int(timestep_value)

        ax.set_title(f"Timestep: {timestep_value}")
        return (scat,)

    num_frames = len(selected_positions) // frame_skip

    ani = FuncAnimation(fig, update, frames=num_frames, interval=interval)

    # Show the animation
    plt.show()


# Example usage:
filename = "code_p3/animation_0.7/int_f_0.7_w_v_0.7_.txt"
first_timestep = 18000
last_timestep = 25000
fps = 30
frame_skip = 10
zoom = 500

animate_particles(
    filename, first_timestep, last_timestep, fps=fps, frame_skip=frame_skip, zoom=zoom
)
