import pyvista as pv
import os

def ply_to_png(ply_path, png_path):
    # read .ply
    mesh = pv.read(ply_path)
    # render & save .png
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(mesh, color= (135/255, 206/255, 235/255))
    plotter.camera_position = [(350, 280, 370), (0, 0, 0), (0, 1, 0)]
    plotter.reset_camera_clipping_range()
    plotter.screenshot(png_path)
    plotter.close()

input_folder = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/ply"
output_folder = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/png"

# 遍历文件夹中的所有.ply文件
for file in os.listdir(input_folder):
    if file.endswith(".ply"):
        ply_file_path = os.path.join(input_folder, file)
        png_file_path = os.path.join(output_folder, os.path.splitext(file)[0] + '.png')
        ply_to_png(ply_file_path, png_file_path)



