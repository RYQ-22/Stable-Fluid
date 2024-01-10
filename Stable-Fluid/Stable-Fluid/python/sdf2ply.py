import numpy as np
from skimage import measure
from plyfile import PlyData, PlyElement
import argparse
import open3d as o3d

def read_sdf(input_file):
    with open(input_file, 'r') as file:
        # read the first row
        line = file.readline().strip().split()
        N1, N2, N3, voxel_size = map(float, line)
        N1 = int(N1)
        N2 = int(N2)
        N3 = int(N3)
        # init np array
        sdf_array = np.zeros((N1, N2, N3))

        # read sdf by rows
        for i in range(N1):
            for j in range(N2):
                for k in range(N3):
                    sdf_value = float(file.readline().strip())
                    sdf_array[i, j, k] = sdf_value

    return sdf_array, voxel_size

def marching_cubes(sdf, voxel_size, threshold=0.0):
    verts, faces, normals, values = measure.marching_cubes(sdf, level=threshold, spacing=(voxel_size, voxel_size, voxel_size))
    return verts, faces

def save_ply(vertices, faces, output_file):
    vertex_element = PlyElement.describe(vertices, 'vertex')
    face_element = PlyElement.describe(faces, 'face')

    ply_data = PlyData([vertex_element, face_element])
    ply_data.write(output_file)
    # print(f"{output_file} is saved successfully.\n")

def main():
    parser = argparse.ArgumentParser(description='Process some arguments.')
    parser.add_argument('-frame', type=int, help='Frame number')
    args = parser.parse_args()
    frame_number = args.frame
    # generate model
    sdf, voxel_size = read_sdf(f"C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/sdf/sdf_{frame_number}.txt")

    

    # sdf to mesh
    vertices, faces = marching_cubes(sdf, voxel_size)
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(vertices)
    mesh.triangles = o3d.utility.Vector3iVector(faces)
    mesh = mesh.filter_smooth_simple(number_of_iterations=2)
    mesh = mesh.filter_smooth_laplacian(number_of_iterations=2)
    mesh = mesh.subdivide_loop(number_of_iterations=2)
    vertices = np.asarray(mesh.vertices)
    faces = np.asarray(mesh.triangles)
    vertices *= 10
    dtype1 = [('x', 'f4'), ('y', 'f4'), ('z', 'f4')]
    structured_vertices = np.array([tuple(vertex) for vertex in vertices], dtype=dtype1)    
    dtype2 = [('vertex_indices', 'i4', (3,))]    
    structured_faces = np.array([(face,) for face in faces], dtype=dtype2)

    # save as .ply file
    save_ply(structured_vertices, structured_faces, f"C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/ply/mesh_{frame_number}.ply")

if __name__ == "__main__":
    main()
    