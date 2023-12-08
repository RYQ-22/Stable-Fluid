from pysdf import SDF
import trimesh
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-input_file', type=str, help='Description for input file')
    parser.add_argument('-output_file', type=str, help='Description for input file')
    parser.add_argument('-N1', type=int, help='Description for N1', default=20)
    parser.add_argument('-N2', type=int, help='Description for N2', default=20)
    parser.add_argument('-N3', type=int, help='Description for N3', default=20)
    parser.add_argument('-size', type=int, help='size', default=20)
    parser.add_argument('-voxel_size', type=float, help='Voxel size', default=0.1)
    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    N1 = args.N1
    N2 = args.N2
    N3 = args.N3
    size = args.size
    voxel_size = args.voxel_size

    # generate np array
    # pos: center of ceil  
    x = np.linspace(-voxel_size*N1/2, voxel_size*N1/2, num=N1)
    y = np.linspace(-voxel_size*N2/2, voxel_size*N2/2, num=N2)
    z = np.linspace(-voxel_size*N3/2, voxel_size*N3/2, num=N3)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    pos = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T

    mesh = trimesh.load(input_file)
    # translate grid    
    vertices = mesh.vertices
    min_xyz = np.amin(vertices, axis=0)
    max_xyz = np.amax(vertices, axis=0)
    delta_x = abs(max_xyz[0] - min_xyz[0])
    mid_x = (max_xyz[0] + min_xyz[0])/2
    mid_y = (max_xyz[1] + min_xyz[1])/2
    mid_z = (max_xyz[2] + min_xyz[2])/2
    offsets = np.array([mid_x, mid_y, mid_z])
    pos *= delta_x / (size * voxel_size)
    pos += offsets
    # compute sdf
    f = SDF(mesh.vertices, mesh.faces)
    sdf = -f(pos)
    
    sdf *= size * voxel_size / delta_x
    np.savetxt(output_file, sdf, header=f'{N1} {N2} {N3} {voxel_size}', comments='')
    # print(f"{output_file} is saved successfully.\n")

if __name__ == "__main__":
    main()