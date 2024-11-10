from moviepy.editor import VideoFileClip

# Load the AVI file
clip = VideoFileClip("wave_plane_animation.avi")

# Convert to GIF and save
clip.write_gif("wave_plane_animation.gif", fps=100)
