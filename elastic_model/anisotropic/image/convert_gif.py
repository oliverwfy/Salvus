from moviepy.editor import VideoFileClip

# Load the AVI file
clip = VideoFileClip("2_layers.mp4")

# Convert to GIF and save
clip.write_gif("2_layers.gif", fps=100)
