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
    parser.add_argument('-size', type=int, help='size')
    parser.add_argument('-scale', type=float, help='scale')
    parser.add_argument('-voxel_size', type=float, help='Voxel size', default=0.1)
    parser.add_argument('-translate_x', type=float, help='offset x', default=0.0)
    parser.add_argument('-translate_y', type=float, help='offset y', default=0.0)
    parser.add_argument('-translate_z', type=float, help='offset z', default=0.0)
    parser.add_argument('-rotate', type=bool, help='whether rotate', default=False)

    args = parser.parse_args()
    size = 0
    scale = 0

    use_size = True
    if args.size is not None and args.scale is not None:
        parser.error("Please provide either -size or -scale, not both.")
    elif args.size is not None:
        use_size = True
        size = args.size
    elif args.scale is not None:
        use_size = False
        scale = args.scale
    else:
        parser.error("Please provide either -size or -scale.")

    input_file = args.input_file
    output_file = args.output_file
    N1 = args.N1
    N2 = args.N2
    N3 = args.N3
    voxel_size = args.voxel_size
    translate_x = args.translate_x
    translate_y = args.translate_y
    translate_z = args.translate_z
    rotate = args.rotate

    if use_size:
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
        if rotate:
            vertices[:, [0, 2]] = vertices[:, [2, 0]]
            vertices[:, 0] = -vertices[:, 0]
            mesh.vertices = vertices
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
    else:
        x = np.linspace(0, voxel_size*N1, num=N1)
        y = np.linspace(0, voxel_size*N2, num=N2)
        z = np.linspace(0, voxel_size*N3, num=N3)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        pos = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T

        mesh = trimesh.load(input_file)
        # translate grid    
        vertices = mesh.vertices
        offsets = np.array([translate_x, translate_y, translate_z])
        pos -= offsets        
        pos /= scale
        # compute sdf
        f = SDF(mesh.vertices, mesh.faces)
        sdf = -f(pos)

        sdf *= scale
        np.savetxt(output_file, sdf, header=f'{N1} {N2} {N3} {voxel_size}', comments='')

if __name__ == "__main__":
    main()